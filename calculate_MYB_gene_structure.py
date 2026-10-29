#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import gzip
import json
from pathlib import Path
from urllib.parse import unquote

METRICS = ['Genelen', 'ExonLen', 'ExonNum', 'CDSLen', 'CDSNum',
           'IntronLen', 'IntronNum', 'AALen']

def text_open(path):
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8-sig') if path.suffix == '.gz' else path.open(encoding='utf-8-sig')

def fasta_records(path):
    records = {}
    identifier = None
    pieces = []
    def finish():
        if identifier is None:
            return
        if identifier in records:
            raise ValueError(f'Duplicate FASTA identifier: {identifier}')
        sequence = ''.join(pieces)
        if not sequence:
            raise ValueError(f'Empty FASTA sequence: {identifier}')
        records[identifier] = sequence
    with text_open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                finish()
                identifier = line[1:].split()[0]
                pieces = []
            else:
                if identifier is None or any(c.isspace() for c in line):
                    raise ValueError('Malformed FASTA sequence')
                pieces.append(line)
    finish()
    return records

def parse_attributes(text):
    return {k: unquote(v) for k,v in
            (x.split('=',1) for x in text.rstrip(';').split(';') if '=' in x)}

def parse_gff(path, prefix='', selected_genes=None):
    genes, parents = {}, {}
    features = []
    origin = Path(path).name
    with text_open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            columns = line.rstrip('\r\n').split('\t')
            if len(columns) != 9:
                raise ValueError(f'{origin}:{lineno}: expected 9 GFF3 columns')
            attrs = parse_attributes(columns[8])
            start, end = int(columns[3]), int(columns[4])
            if start > end or start < 1:
                raise ValueError(f'{origin}:{lineno}: invalid coordinates')
            kind = columns[2]
            if kind == 'gene':
                key = prefix + attrs['ID']
                if key in genes:
                    raise ValueError(f'Duplicate GFF gene identifier: {key}')
                genes[key] = {'start':start, 'end':end, 'seqid':columns[0],
                              'strand':columns[6], 'exons':[], 'cds':[],
                              'source':origin, 'transcripts':[]}
            elif kind in ('mRNA', 'transcript'):
                key = prefix + attrs['ID']
                if key in parents:
                    raise ValueError(f'Duplicate transcript identifier: {key}')
                parents[key] = [prefix+x for x in attrs.get('Parent','').split(',') if x]
            elif kind in ('exon', 'CDS'):
                parent_ids = [prefix+x for x in attrs.get('Parent','').split(',') if x]
                features.append((kind, start, end, columns[0], columns[6], parent_ids, lineno))
    def resolve(key, seen=None):
        if key in genes:
            return {key}
        seen = set() if seen is None else set(seen)
        if key in seen:
            raise ValueError(f'GFF Parent cycle at {key}')
        seen.add(key)
        result = set()
        for parent in parents.get(key, []):
            result.update(resolve(parent, seen))
        return result
    for transcript in parents:
        for gene in resolve(transcript):
            genes[gene]['transcripts'].append(transcript)
    unresolved = []
    for kind, start, end, seqid, strand, parent_ids, lineno in features:
        gene_ids = set()
        for parent in parent_ids:
            gene_ids.update(resolve(parent))
        if not gene_ids and selected_genes is None:
            unresolved.append({'source':origin, 'line':lineno, 'feature':kind,
                               'parents':','.join(parent_ids)})
        for key in gene_ids:
            if selected_genes is not None and key not in selected_genes:
                continue
            gene = genes[key]
            if seqid != gene['seqid'] or strand not in (gene['strand'], '.'):
                raise ValueError(f'{origin}:{lineno}: feature/gene chromosome or strand mismatch')
            if start < gene['start'] or end > gene['end']:
                raise ValueError(f'{origin}:{lineno}: feature lies outside its gene')
            gene['exons' if kind == 'exon' else 'cds'].append((start,end))
    if selected_genes is not None:
        genes={key:value for key,value in genes.items() if key in selected_genes}
    return genes, unresolved

def calculate(gene, protein):
    if len(gene['transcripts']) != 1:
        raise ValueError(f'Expected one retained transcript, found {len(gene["transcripts"])}')
    cds = sorted(gene['cds'])
    if not cds:
        raise ValueError('No CDS features')
    exons = sorted(gene['exons'] or cds)
    for name, intervals in [('exon',exons),('CDS',cds)]:
        if len(intervals) != len(set(intervals)):
            raise ValueError(f'Duplicate {name} intervals')
        if any(right[0] <= left[1] for left,right in zip(intervals,intervals[1:])):
            raise ValueError(f'Overlapping {name} intervals')
    gene_length = gene['end']-gene['start']+1
    exon_length = sum(end-start+1 for start,end in exons)
    gap_length = sum(right[0]-left[1]-1 for left,right in zip(exons,exons[1:]))
    return {'Genelen':gene_length, 'ExonLen':exon_length, 'ExonNum':len(exons),
            'CDSLen':sum(end-start+1 for start,end in cds), 'CDSNum':len(cds),
            'IntronLen':gene_length-exon_length, 'IntronNum':max(0,len(exons)-1),
            'AALen':len(protein), 'IntronGapLen':gap_length,
            'TerminalNonExonicLen':gene_length-exon_length-gap_length,
            'ExonBasis':'annotated_exon' if gene['exons'] else 'CDS_proxy_no_exon_feature',
            'AnnotationSource':gene['source']}

def read_tsv(path):
    with text_open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write_tsv(path, fieldnames, rows):
    with Path(path).open('w',encoding='utf-8',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=fieldnames,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def run(config_path, output_dir, strict=False):
    config_path = Path(config_path).resolve()
    folder = config_path.parent
    config = json.loads(config_path.read_text(encoding='utf-8'))
    return run_calculation(config,folder,output_dir,strict)

def run_calculation(config, folder, output_dir, strict=False):
    def source_path(value):
        p = Path(value)
        return p if p.is_absolute() else folder/p
    selected=set(config['selected_genes']) if config.get('selected_genes') else None
    genes, unresolved = parse_gff(source_path(config['gff']),config.get('gff_prefix',''),selected)
    overlay_log = []
    for item in config.get('annotation_overlays',[]):
        overlay, warnings = parse_gff(source_path(item['path']),item.get('prefix',''))
        unresolved.extend(warnings)
        for key, gene in overlay.items():
            if selected is not None and key not in selected:
                continue
            if key not in genes:
                raise ValueError(f'Overlay gene absent from main annotation: {key}')
            old = genes[key]
            if (gene['start'],gene['end'],gene['strand']) != (old['start'],old['end'],old['strand']):
                raise ValueError(f'Overlay changes the genomic bounds of {key}')
            genes[key] = gene
            overlay_log.append({'GFFGeneID':key,'source':item['path'],'reason':item['reason']})
    proteins = fasta_records(source_path(config['protein_fasta']))
    if config.get('gene_map'):
        mapping = read_tsv(source_path(config['gene_map']))
    else:
        mapping=[]
        with text_open(source_path(config['protein_fasta'])) as handle:
            for line in handle:
                if not line.startswith('>'): continue
                parts=line[1:].strip().split()
                attrs=dict(p.split('=',1) for p in parts[1:] if '=' in p)
                mapping.append({'GeneID':parts[0], 'GeneName':attrs.get('GeneName',parts[0]),
                                'GFFGeneID':attrs.get('source_id',parts[0])})
    if len({r['GeneID'] for r in mapping}) != len(mapping):
        raise ValueError('Nonunique output GeneID in map')
    if len({r['GFFGeneID'] for r in mapping}) != len(mapping):
        raise ValueError('Nonunique source gene in map')
    if selected is not None:
        mapping=[r for r in mapping if r['GFFGeneID'] in selected]
        missing=selected-{r['GFFGeneID'] for r in mapping}
        if missing:
            raise ValueError('Requested genes absent from protein/mapping input: '+','.join(sorted(missing)))
    if not mapping:
        raise ValueError('No genes selected')
    rows, errors = [], []
    for item in mapping:
        key = item['GFFGeneID']
        protein_key = key if key in proteins else item['GeneID']
        try:
            if key not in genes or protein_key not in proteins:
                raise ValueError('Gene absent from GFF or protein FASTA')
            values = calculate(genes[key],proteins[protein_key])
        except ValueError as error:
            if not config.get('keep_going'):
                raise ValueError(f'{key}: {error}') from error
            errors.append({'GeneID':item['GeneID'],'GFFGeneID':key,'Error':str(error)})
            continue
        rows.append({'GeneID':item['GeneID'],'GeneName':item['GeneName'],
                     'GFFGeneID':key,**values})
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True,exist_ok=True)
    fields=['GeneID','GeneName','GFFGeneID']+METRICS+['IntronGapLen','TerminalNonExonicLen','ExonBasis','AnnotationSource']
    write_tsv(output_dir/'recomputed_gene_structure.tsv',fields,rows)
    write_tsv(output_dir/'input_errors.tsv',['GeneID','GFFGeneID','Error'],errors)
    write_tsv(output_dir/'unresolved_main_GFF_features.tsv',
              ['source','line','feature','parents'],unresolved)
    discrepancies = []
    reference = ({r['GeneID']:r for r in read_tsv(source_path(config['reference_metrics']))}
                 if config.get('reference_metrics') else {})
    for row in rows:
        if not reference: break
        expected = reference[row['GeneID']]
        for metric in METRICS:
            if row[metric] != int(expected[metric]):
                discrepancies.append({'GeneID':row['GeneID'],'GeneName':row['GeneName'],
                                      'metric':metric,'reference':int(expected[metric]),
                                      'recomputed':row[metric]})
    write_tsv(output_dir/'reference_discrepancies.tsv',
              ['GeneID','GeneName','metric','reference','recomputed'],discrepancies)
    report = {'status':('matches_reference' if not discrepancies else 'reference_discrepancies_present')
                         if reference else 'calculated_without_reference_comparison',
              'genes_selected':len(mapping),'genes_recomputed':len(rows),'genes_with_input_errors':len(errors),
              'reference_metric_cells':len(rows)*len(METRICS) if reference else 0,
              'matching_genes':len(rows)-len({r['GeneID'] for r in discrepancies}) if reference else None,
              'mismatching_genes':len({r['GeneID'] for r in discrepancies}),
              'mismatching_metric_cells':len(discrepancies),
              'annotation_overlays':overlay_log,'unresolved_main_features':len(unresolved),
              'exon_CDS_proxy_genes':sum(r['ExonBasis'].startswith('CDS_proxy') for r in rows),
              'genes_with_terminal_non_exonic_bases':sum(r['TerminalNonExonicLen']>0 for r in rows)}
    if errors:
        report['status']='input_errors_present'
    (output_dir/'verification_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    return 2 if errors or (strict and discrepancies) else 0

def main():
    parser = argparse.ArgumentParser(description='Calculate gene-structure metrics from GFF3 and protein FASTA.')
    parser.add_argument('--config',type=Path,help='Optional audit configuration with exact input mapping and annotation overlays')
    parser.add_argument('--gff',type=Path,help='GFF3 containing the retained gene models')
    parser.add_argument('--proteins',type=Path,help='Protein FASTA (source_id= header fields supported)')
    parser.add_argument('--gene',action='append',help='Calculate only this exact GFF gene ID; repeat for multiple genes')
    parser.add_argument('--gff-prefix',default='',help='Prefix to add to GFF IDs before exact matching')
    parser.add_argument('--keep-going',action='store_true',help='Report per-gene input errors and continue with other genes')
    parser.add_argument('--gene-map',type=Path,help='Optional TSV: GeneID, GeneName, GFFGeneID')
    parser.add_argument('--reference',type=Path,help='Optional TSV: GeneID and the eight metric columns; validation only')
    parser.add_argument('--out',type=Path,default=Path('reproduced_results'))
    parser.add_argument('--strict',action='store_true',help='Exit 2 if any metric differs from the supplied reference table')
    args = parser.parse_args()
    if args.config:
        if args.gff or args.proteins or args.gene_map or args.reference or args.gene or args.gff_prefix or args.keep_going:
            parser.error('--config cannot be combined with direct input arguments')
        code=run(args.config,args.out,args.strict)
    else:
        if not args.gff or not args.proteins:
            parser.error('provide --gff and --proteins, or an explicit --config')
        config={'gff':str(args.gff.resolve()),'protein_fasta':str(args.proteins.resolve())}
        if args.gene: config['selected_genes']=args.gene
        if args.gff_prefix: config['gff_prefix']=args.gff_prefix
        if args.keep_going: config['keep_going']=True
        if args.gene_map: config['gene_map']=str(args.gene_map.resolve())
        if args.reference: config['reference_metrics']=str(args.reference.resolve())
        code=run_calculation(config,Path.cwd(),args.out,args.strict)
    raise SystemExit(code)

if __name__ == '__main__':
    main()

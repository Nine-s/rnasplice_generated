nextflow.enable.dsl = 2


include { FASTQC  } from '../../../data/ninon/description_prototype/modules_nextflow/fastqc.nf'
include { TRIMGALORE  } from '../../../data/ninon/description_prototype/modules_nextflow/trimgalore.nf'
include { SALMON_GENOMEGENERATE  } from '../../../data/ninon/description_prototype/modules_nextflow/salmon_genome_generate.nf'
include { SALMON_QUANT  } from '../../../data/ninon/description_prototype/modules_nextflow/salmon.nf'
include { STAR_GENOMEGENERATE  } from '../../../data/ninon/description_prototype/modules_nextflow/star_genome_generate.nf'
include { STAR_ALIGN  } from '../../../data/ninon/description_prototype/modules_nextflow/star_align.nf'
include { SAMTOOLS  } from '../../../data/ninon/description_prototype/modules_nextflow/samtools.nf'
include { CUSTOM_GETCHROMSIZES  } from '../../../data/ninon/description_prototype/modules_nextflow/getchromsizes.nf'
include { BEDTOOLS_GENOMECOV  } from '../../../data/ninon/description_prototype/modules_nextflow/bedtoolsgenomecov.nf'
include { BEDCLIP as BEDCLIP_FORWARD; BEDCLIP as BEDCLIP_REVERSE } from '../../../data/ninon/description_prototype/modules_nextflow/bedclip.nf'
include { BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD; BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE } from '../../../data/ninon/description_prototype/modules_nextflow/bedgraphtobigwig.nf'
include { DEXSEQ_ANNOTATION  } from '../../../data/ninon/description_prototype/modules_nextflow/dexseq_annotation.nf'
include { DEXSEQ_COUNT  } from '../../../data/ninon/description_prototype/modules_nextflow/dexseq_count.nf'
include { MERGE_RESULTS_DEXSEQ  } from '../../../data/ninon/description_prototype/modules_nextflow/merge_results_dexseq.nf'
include { DEXSEQ_EXON  } from '../../../data/ninon/description_prototype/modules_nextflow/dexseq_exon.nf'
include { GFFREAD_TX2GENE  } from '../../../data/ninon/description_prototype/modules_nextflow/gffread_tx2gene.nf'
include { MERGE_RESULTS_SALMON  } from '../../../data/ninon/description_prototype/modules_nextflow/merge_results.nf'
include { TXIMPORT  } from '../../../data/ninon/description_prototype/modules_nextflow/tximport.nf'
include { DRIMSEQ_FILTER  } from '../../../data/ninon/description_prototype/modules_nextflow/drimseq_filter.nf'
include { DEXSEQ_DTU ; MULTIQC  } from '../../../data/ninon/description_prototype/modules_nextflow/dexseq_dtu.nf'

workflow{
        read_pairs_ch = Channel
            .fromPath( params.csv_input )
            .splitCsv(header: true, sep: ',')
            .map {row -> tuple(row.sample, [row.path_r1, row.path_r2])}
            .view()
        
FASTQC(params.samples)
SALMON_GENOMEGENERATE(params.genome, params.transcripts_fasta)
TRIMGALORE(params.samples)
DEXSEQ_ANNOTATION(params.annotation_gtf)
CUSTOM_GETCHROMSIZES(params.genome)
GFFREAD_TX2GENE(params.annotation_gtf)
STAR_GENOMEGENERATE(params.genome, params.annotation_gtf)
SALMON_QUANT(TRIMGALORE.out.preprocessed_reads, SALMON_GENOMEGENERATE.out.index)
MERGE_RESULTS_SALMON(SALMON_QUANT.out.transcripts.collect())
STAR_ALIGN(TRIMGALORE.out.preprocessed_reads, STAR_GENOMEGENERATE.out.index, params.annotation_gtf)
SAMTOOLS(STAR_ALIGN.out.sam)
TXIMPORT(MERGE_RESULTS_SALMON.out.gathered_bam, GFFREAD_TX2GENE.out.tx2gene)
DEXSEQ_COUNT(SAMTOOLS.out.bam, DEXSEQ_ANNOTATION.out.gff, params.alignment_quality)
MULTIQC(SALMON_QUANT.out.json_info.collect(), TRIMGALORE.out.log.collect(), STAR_ALIGN.out.log_final.collect(), FASTQC.out.zip.collect())
BEDTOOLS_GENOMECOV(SAMTOOLS.out.bam)
BEDCLIP_FORWARD(BEDTOOLS_GENOMECOV.out.bedgraph_forward, CUSTOM_GETCHROMSIZES.out.sizes)
DRIMSEQ_FILTER(TXIMPORT.out.txi_dtu, TXIMPORT.out.tximport_tx2gene, params.csv_input, params.min_samps_gene_expr, params.min_samps_feature_expr, params.min_samps_feature_prop, params.min_feature_expr, params.min_feature_prop, params.min_gene_expr)
MERGE_RESULTS_DEXSEQ(DEXSEQ_COUNT.out.dexseq_clean_txt.collect())
DEXSEQ_DTU(DRIMSEQ_FILTER.out.drimseq_samples_tsv, DRIMSEQ_FILTER.out.drimseq_counts_tsv, params.csv_contrastsheet, params.n_dexseq_plot)
DEXSEQ_EXON(MERGE_RESULTS_DEXSEQ.out.clean_counts, DEXSEQ_ANNOTATION.out.gff, params.csv_input, params.csv_contrastsheet, params.n_dexseq_plot)
BEDCLIP_REVERSE(BEDTOOLS_GENOMECOV.out.bedgraph_reverse, CUSTOM_GETCHROMSIZES.out.sizes)
BEDGRAPH_TO_BIGWIG_FORWARD(BEDCLIP_FORWARD.out.bedgraph, CUSTOM_GETCHROMSIZES.out.sizes)
BEDGRAPH_TO_BIGWIG_REVERSE(BEDCLIP_REVERSE.out.bedgraph, CUSTOM_GETCHROMSIZES.out.sizes)

}
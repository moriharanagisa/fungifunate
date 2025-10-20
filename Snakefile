import os
configfile: "config/config.yaml"

reads    = config.get("reads")
annot    = config.get("annot")
meta     = config.get("meta")
outdir   = config.get("outdir")

# 入力推定
longest_pep = os.path.join(annot, "transcripts.fasta.transdecoder_dir/longest_orfs.pep")
longest_cds = os.path.join(annot, "transcripts.fasta.transdecoder_dir/longest_orfs.cds")
ggsearch_txt = expand(os.path.join(annot, "ggsearch-{db}.txt"),
                      db=["shiitake-human","shiitake-mouse","shiitake-S.cerevisiae","shiitake-uniprot","shiitake-fungidb"])
interpro_tsv = os.path.join(annot, "interpro_longest_orfs.tsv")
uniprot_go   = os.path.join(annot, "uniprot_go_annotations.tsv")

SAMPLES = []
for f in sorted(os.listdir(reads)):
    if f.endswith("_1.fastq"):
        srr = f.split("_1.fastq")[0]
        if os.path.exists(os.path.join(reads, f"{srr}_2.fastq")):
            SAMPLES.append(srr)

rule all:
    input:
        os.path.join(outdir, "annotation/ggsearch_annotbl.txt"),
        os.path.join(outdir, "annotation/ggsearch-interpro-GO_annotbl.txt"),
        expand(os.path.join(outdir, "salmon", "salmon-{s}/quant.sf"), s=SAMPLES),
        os.path.join(outdir, "tximport", "tpm.tsv"),
        os.path.join(outdir, "de", "deseq2-primordia-vs-mycelia.txt"),
        os.path.join(outdir, "pca", "PCA_plot.png"),
        os.path.join(outdir, "topgo", "topGO_dotplot_weight01.png")

rule get_longest_orf_pid:
    input: longest_pep
    output: os.path.join(outdir, "annotation/longest_orf-pid.txt")
    shell: r"""
      mkdir -p {outdir}/annotation
      perl -nle 'print $1 if(/^>(\S+)/)' {input} > {output}
    """

rule parse_ggsearch_each:
    input: txt=ggsearch_txt
    output: parsed=expand(os.path.join(outdir,"annotation","{b}-parsing.txt"), b=[os.path.basename(x).replace(".txt","") for x in ggsearch_txt])
    shell: r"""
      mkdir -p {outdir}/annotation
      for file in {annot}/ggsearch-*.txt; do
        perl {workflow.basedir}/bin/15parseggsearch.pl < "$file" > {outdir}/annotation/$(basename "$file" .txt)-parsing.txt
      done
    """

rule make_annotation_table:
    input:
        pid = os.path.join(outdir,"annotation/longest_orf-pid.txt"),
        parsed = lambda wc: expand(os.path.join(outdir,"annotation","ggsearch-shiitake-{db}-parsing.txt"),
                                   db=["human","mouse","S.cerevisiae","uniprot","fungidb"])
    output: os.path.join(outdir,"annotation/ggsearch_annotbl.txt")
    shell: r"""
      cat {input.pid} \
      | perl {workflow.basedir}/bin/15mkannotbl-for-fungi.pl {outdir}/annotation/ggsearch-shiitake-human-parsing.txt \
      | perl {workflow.basedir}/bin/15mkannotbl.pl {outdir}/annotation/ggsearch-shiitake-mouse-parsing.txt \
      | perl {workflow.basedir}/bin/15mkannotbl.pl {outdir}/annotation/ggsearch-shiitake-S.cerevisiae-parsing.txt \
      | perl {workflow.basedir}/bin/15mkannotbl.pl {outdir}/annotation/ggsearch-shiitake-uniprot-parsing.txt \
      | perl {workflow.basedir}/bin/15mkannotbl.pl {outdir}/annotation/ggsearch-shiitake-fungidb-parsing.txt \
      > {output}
    """

rule join_interpro_and_go:
    input:
        annot = os.path.join(outdir,"annotation/ggsearch_annotbl.txt"),
        interpro = interpro_tsv,
        uniprot_go = uniprot_go
    output:
        os.path.join(outdir,"annotation/ggsearch-interpro-GO_annotbl.txt")
    shell: r"""
      Rscript {workflow.basedir}/bin/ggsearch_join.R \
        --annot {input.annot} \
        --interpro {input.interpro} \
        --uniprot {input.uniprot_go} \
        --out {output}
    """

rule salmon_index:
    input: longest_cds
    output: directory(os.path.join(outdir,"salmon","index"))
    params: k="--keepDuplicates"
    threads: 8
    shell: r"""
      mkdir -p {outdir}/salmon
      salmon index -p {threads} {params.k} -t {input} -i {output}
    """

rule salmon_quant:
    input:
        idx = os.path.join(outdir,"salmon/index"),
        r1  = lambda wc: os.path.join(reads, f"{wc.s}_1.fastq"),
        r2  = lambda wc: os.path.join(reads, f"{wc.s}_2.fastq")
    output:
        os.path.join(outdir,"salmon","salmon-{s}","quant.sf")
    threads: 8
    shell: r"""
      salmon quant -i {input.idx} -l A \
        -1 {input.r1} -2 {input.r2} \
        -p {threads} -o {outdir}/salmon/salmon-{wildcards.s}
    """

rule tximport_deseq_topgo_all:
    input:
        meta = meta,
        annot = os.path.join(outdir,"annotation/ggsearch-interpro-GO_annotbl.txt"),
        quants = expand(os.path.join(outdir,"salmon","salmon-{s}","quant.sf"), s=SAMPLES)
    output:
        tpm = os.path.join(outdir,"tximport/tpm.tsv"),
        deg1 = os.path.join(outdir,"de/deseq2-primordia-vs-mycelia.txt"),
        deg2 = os.path.join(outdir,"de/deseq2-fruiting_body-vs-primordia.txt"),
        deg3 = os.path.join(outdir,"de/deseq2-mycelia-vs-fruiting_body.txt"),
        pca = os.path.join(outdir,"pca/PCA_plot.png"),
        topgo = os.path.join(outdir,"topgo/topGO_dotplot_weight01.png")
    threads: 8
    shell: r"""
      Rscript {workflow.basedir}/bin/analysis.R \
        --meta {input.meta} \
        --annot {input.annot} \
        --salmon_dir {outdir}/salmon \
        --out {outdir}
    """

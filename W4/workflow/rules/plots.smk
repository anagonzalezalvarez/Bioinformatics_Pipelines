#Generate UMAP and tSNE plots 

rule plots:
    input:
        get_data_path
    output:
        umap_plot = os.path.join(result_path,'{analysis}','plots','UMAP_plot.png'),
        tsne_plot = os.path.join(result_path,'{analysis}','plots','tSNE_plot.png'),                                      
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","{analysis}_plots.log"),
    params:
        partition=config.get("partition"),
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
    script:
        "../scripts/plots.R"


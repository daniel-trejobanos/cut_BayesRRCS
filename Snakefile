rule maf_sd:
    input:
        sd = "data/ukb_imp_v3_UKB_EST_oct19_unrelated_full_stats.rds",
        maf = "data/ukb_imp_v3_UKB_EST_oct19_unrelated.frq"
    output:
        maf_sd = "data/maf_sd.rds"
    shell:
        "Rscript R/generate_maf_sds.R --maf {input.maf} --sd {input.sd} --out {output}"

rule snp_groups:
    input:
        group_id="data/ukb_imp_v3_UKB_EST_oct19_unrelated_78groups_13annot_3maf_2ld.group",
        snp_id="data/snplist.txt"
    output:
        "data/snp_groups.rds"
    shell:
        "Rscript R/create_groups.R --snp_id {input.snp_id} --groups {input.group_id} --out {output}"
            
rule test_unscale:
    input:
        chain = "data/groups78_mix4_cpus4_tasks8_nodes10_mpisync10_height_1_unrelated_148covadjusted_w35520NA.betMap",
        mafsd = rules.maf_sd.output
    output:
        "test/unscaled_test.rds"
    shell:
        "Rscript R/unscale_betas.R --chain {input.chain} --start 500 --end 501 --mafsd {input.mafsd} --out {output}"


```mermaid

flowchart TD
    R["raw reads\nvrouw_maria_2026_segments/"] --> S1

    subgraph S1["① bwa-vrouw.sh"]
        direction TB
        A1["SAMPLE_1.sai + _2.sai\n~1GB"]
        A2["bwa-out/SAMPLE_sorted.bam\n12GB"]
    end
    S1 --> S2

    subgraph S2["② add_readgroups.sh"]
        B1["bwa-out/SAMPLE_sorted_RG.bam\n13GB"]
    end
    S2 --> S3

    subgraph S3["③ mapDamage_vrouw_maria.sh"]
        direction TB
        C1["bwa-out/SAMPLE_dedup.bam\n13GB"]
        C2["mapDamage-out/SAMPLE_mapDamage/\nQC plots, few MB"]
    end
    S3 --> S4

    subgraph S4["④ mapDamage_rescale.sh"]
        D1["mapDamage-out/.../SAMPLE_dedup.rescaled.bam\n12GB"]
    end
    S4 --> S5

    subgraph S5["⑤ gatk_hc.sh"]
        E1["gatk-out/SAMPLE.g.vcf.gz\n50KB"]
    end
    S5 --> S6

    subgraph S6["⑥ genotype_gvcf.sh"]
        F1["gatk-out/SAMPLE.vcf.gz\n11KB"]
    end
    S6 --> S7

    subgraph S7["⑦ gatk-filter.sh"]
        G1["gatk-out/SAMPLE_filtered.vcf.gz\n11KB"]
    end
    S7 --> VIZ["visualization/\nextract_snp_viz_data.sh"]

    classDef used fill:#1a7f1a,color:#fff,stroke:#0d4d0d,stroke-width:2px
    classDef dead fill:#d8d8d0,color:#333,stroke:#999
    class A1,A2,B1,C1 dead
    class C2,D1,E1,F1,G1 used

    ```


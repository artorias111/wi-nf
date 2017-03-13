#!/usr/bin/env nextflow
tmpdir = config.tmpdir
reference = config.reference
annotation_reference = config.annotation_reference
alignment_cores = config.alignment_cores
variant_cores = config.variant_cores
genome = config.genome
analysis_dir = config.analysis_dir
SM_alignments_dir = config.SM_alignments_dir
beagle_location = config.beagle_location

// Define contigs here!
contig_list = ["I", "II", "III", "IV", "V", "X", "MtDNA"]
contigs = Channel.from(contig_list)

/*
    Filtering configuration
*/

min_depth=10
qual=30
mq=40
dv_dp=0.5

println "Running Wild Isolate Pipeline!!!"
println "Using Reference: ${genome}" 

// Construct strain and isotype lists
import groovy.json.JsonSlurper

if (test == 'true') {
    strain_json = "isotype_set_test.json"
    min_depth=1
    qual=1
    mq=1
    dv_dp=0
} else {
    strain_json = "isotype_set.json"
}

def strain_set = []

def strainFile = new File(strain_json)
def strainJSON = new JsonSlurper().parseText(strainFile.text)
strain_set_file = Channel.fromPath(strain_json)

strainJSON.each { SM, RG ->
    RG.each { k, v ->
        strain_set << [SM, k, v[0], v[1], v[2]]
    }
}


process setup_dirs {

    executor 'local'

    input:
        file strain_set_file

    """
        mkdir -p ${analysis_dir}
        cp ${strain_set_file} ${analysis_dir}/isotype_set.json
    """
}

/*
    Alignment
*/
process perform_alignment {

    cpus alignment_cores

    tag { fq_pair_id }

    input:
        set SM, RG, fq1, fq2, fq_pair_id from strain_set
    output:
        set val(SM), file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into fq_bam_set

    
    """
        bwa mem -t ${alignment_cores} -R '${RG}' ${reference} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${alignment_cores} --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${alignment_cores} --show-progress --tmpdir=${tmpdir} --out=${fq_pair_id}.bam /dev/stdin
        sambamba index --nthreads=${alignment_cores} ${fq_pair_id}.bam
    """
}

/* 
  Merge - Generate SM Bam
*/

process merge_bam {

    cpus alignment_cores

    publishDir SM_alignments_dir + "/WI/isotype", mode: 'copy', pattern: '*.bam*'

    tag { SM }

    input:
        set SM, bam, index from fq_bam_set.groupTuple()

    output:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") into SM_bam_set 
        file("${SM}.duplicates.txt") into duplicates_file
        
    """
    count=`echo ${bam.join(" ")} | tr ' ' '\\n' | wc -l`

    if [ "\${count}" -eq "1" ]; then
        ln -s ${bam.join(" ")} ${SM}.merged.bam
        ln -s ${bam.join(" ")}.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${alignment_cores} --show-progress ${SM}.merged.bam ${bam.join(" ")}
        sambamba index --nthreads=${alignment_cores} ${SM}.merged.bam
    fi

    picard MarkDuplicates I=${SM}.merged.bam O=${SM}.bam M=${SM}.duplicates.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${alignment_cores} ${SM}.bam
    """
}

SM_bam_set.into { 
                  bam_idxstats; 
                  bam_stats;
                  bam_coverage;
                  bam_snp_individual;
                  bam_snp_union;
                  bam_telseq;
}

process SM_idx_stats {
    
    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_idxstats
    output:
        file bam_idxstats into bam_idxstats_set

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > bam_idxstats
    """
}

process SM_combine_idx_stats {

    publishDir analysis_dir + "/isotype", mode: 'copy'

    input:
        val bam_idxstats from bam_idxstats_set.toSortedList()

    output:
        file("SM_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > SM_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> SM_bam_idxstats.tsv
    """
}

/*
    SM bam stats
*/

process SM_bam_stats {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_stats

    output:
        file 'bam_stat' into SM_bam_stat_files

    """
        cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_SM_bam_stats {

    publishDir analysis_dir + "/SM", mode: 'copy'

    input:
        val stat_files from SM_bam_stat_files.toSortedList()

    output:
        file("SM_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > SM_bam_stats.tsv
        cat ${stat_files.join(" ")} >> SM_bam_stats.tsv
    """
}

process format_duplicates {

    publishDir analysis_dir + "/duplicates", mode: 'copy'

    input:
        val duplicates_set from duplicates_file.toSortedList()

    output:
        file("bam_duplicates.tsv")


    """
        echo -e 'filename\\tlibrary\\tunpaired_reads_examined\\tread_pairs_examined\\tsecondary_or_supplementary_rds\\tunmapped_reads\\tunpaired_read_duplicates\\tread_pair_duplicates\\tread_pair_optical_duplicates\\tpercent_duplication\\testimated_library_size' > bam_duplicates.tsv
        for i in ${duplicates_set.join(" ")}; do
            f=\$(basename \${i})
            cat \${i} | awk -v f=\${f/.duplicates.txt/} 'NR >= 8 && \$0 !~ "##.*" && \$0 != ""  { print f "\\t" \$0 } NR >= 8 && \$0 ~ "##.*" { exit }'  >> bam_duplicates.tsv
        done;
    """
}

/*
    Coverage Bam
*/
process coverage_SM {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_coverage

    output:
        val SM into SM_coverage_sample
        file("${SM}.coverage.tsv") into SM_coverage


    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}

process coverage_SM_merge {

    publishDir analysis_dir + "/isotype", mode: 'copy'

    input:
        val sm_set from SM_coverage.toSortedList()

    output:
        file("SM_coverage.full.tsv")
        file("SM_coverage.tsv") into SM_coverage_merged

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> SM_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > SM_coverage.tsv
    """
}

/*
    telseq
*/

process call_telseq {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_telseq
    output:
        file("telseq_out.txt") into telseq_results

    """
        telseq -z TTAGGC -H ${SM}.bam > telseq_out.txt
    """
}

process combine_telseq {

    publishDir analysis_dir + "/SM"

    input:
        file("telseq?.txt") from telseq_results.toSortedList()

    output:
        file("telseq.tsv")

    """
        telseq -h > telseq.tsv
        cat telseq*.txt >> telseq.tsv
    """
}

/*
    Call Variants
*/

process call_variants_individual {

    tag { SM }

    cpus variant_cores

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_snp_individual

    output:
        file("${SM}.individual.sites.tsv") into individual_sites

    """
    # Perform individual-level calling
    contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
    echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${variant_cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,SP --fasta-ref ${reference} ${SM}.bam | bcftools call --skip-variants indels --variants-only --multiallelic-caller -O z  -  > ${SM}.{}.individual.vcf.gz"
    order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".individual.vcf.gz" }'`
    
    # Output variant sites
    bcftools concat \${order} -O v | vk geno het-polarization - | bcftools view -O z > ${SM}.individual.vcf.gz
    bcftools index ${SM}.individual.vcf.gz
    rm \${order}

    bcftools view -M 2 -m 2 -O v ${SM}.individual.vcf.gz | \\
    bcftools filter --include 'DP > 3' | \\
    egrep '(^#|1/1)' | \\
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' > ${SM}.individual.sites.tsv

    """
}

process merge_variant_list {

    publishDir analysis_dir + "/sitelist", mode: 'copy'
    
    input:
        val sites from individual_sites.toSortedList()

    output:
        set file("sitelist.tsv.gz"), file("sitelist.tsv.gz.tbi") into gz_sitelist
        set file("sitelist.tsv.gz"), file("sitelist.tsv.gz.tbi") into fq_gz_sitelist
        file("sitelist.tsv") into sitelist
        file("sitelist.count.txt")


    """
        echo ${sites}
        cat ${sites.join(" ")} | sort -k1,1 -k2,2n | uniq > sitelist.tsv
        cat sitelist.tsv | wc -l > sitelist.count.txt
        bgzip sitelist.tsv -c > sitelist.tsv.gz && tabix -s1 -b2 -e2 sitelist.tsv.gz
    """
}

union_vcf_set = bam_snp_union.spread(gz_sitelist)


process call_variants_union {

    tag { SM }

    cpus variant_cores

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi') from union_vcf_set

    output:
        val SM into union_vcf_SM
        file("${SM}.union.vcf.gz") into union_vcf_list

    """
        contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
        echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${variant_cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,INFO/AD,SP --fasta-ref ${reference} ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
        order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".union.vcf.gz" }'`
        # "
        # Concatenate and filter
        bcftools concat \${order} -O v | \\
        vk geno het-polarization - | \\
        bcftools filter -O u --threads ${variant_cores} --mode + --soft-filter quality --include "QUAL >= ${qual} || FORMAT/GT == '0/0'" |  \\
        bcftools filter -O u --threads ${variant_cores} --mode + --soft-filter min_depth --include "FORMAT/DP > ${min_depth}" | \\
        bcftools filter -O u --threads ${variant_cores} --mode + --soft-filter mapping_quality --include "INFO/MQ > ${mq}" | \\
        bcftools filter -O u --threads ${variant_cores} --mode + --soft-filter dv_dp --include "(FORMAT/AD[1])/(FORMAT/DP) >= ${dv_dp} || FORMAT/GT == '0/0'" | \\
        bcftools filter --mode + --soft-filter het --exclude 'AC==1' | \\
        vk geno transfer-filter - | \\
        bcftools view -O z > ${SM}.union.vcf.gz
        bcftools index ${SM}.union.vcf.gz
        rm \${order}
    """
}

process generate_union_vcf_list {

    cpus 1 

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
       val vcf_set from union_vcf_list.toSortedList()

    output:
       file("union_vcfs.txt") into union_vcfs

    """
        echo ${vcf_set.join(" ")} | tr ' ' '\\n' > union_vcfs.txt
    """
}

union_vcfs_in = union_vcfs.spread(contigs)

process merge_union_vcf_chromosome {

    cpus alignment_cores

    tag { chrom }

    input:
        set file(union_vcfs:"union_vcfs.txt"), val(chrom) from union_vcfs_in

    output:
        val(chrom) into contigs_list_in
        file("${chrom}.merged.raw.vcf.gz") into raw_vcf

    """
        bcftools merge --threads ${alignment_cores} --regions ${chrom} -O z -m all --file-list ${union_vcfs} > ${chrom}.merged.raw.vcf.gz
        bcftools index ${chrom}.merged.raw.vcf.gz
    """
}

// Generate a list of ordered files.
contig_raw_vcf = contig_list*.concat(".merged.raw.vcf.gz")

process concatenate_union_vcf {

    echo true

    //publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        //val chrom from contigs_list_in
        val merge_vcf from raw_vcf.toSortedList()

    output:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") into raw_vcf_concatenated

    """
        for i in ${merge_vcf.join(" ")}; do
            ln  -s \${i} `basename \${i}`;
        done;
        chrom_set="";
        bcftools concat -O z ${contig_raw_vcf.join(" ")}  > merged.raw.vcf.gz
        bcftools index merged.raw.vcf.gz
    """
}

process filter_merged_vcf {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") from raw_vcf_concatenated

    output:
        set file("filtered.vcf.gz"), file("filtered.vcf.gz.csi") into filtered_vcf
        set val('filtered'), file("filtered.vcf.gz"), file("filtered.vcf.gz.csi") into filtered_vcf_stat

    """
        bcftools view merged.raw.vcf.gz | \\
        vk filter MISSING --max=0.90 --soft-filter="high_missing" --mode=x - | \
        vk filter HET --max=0.10 --soft-filter="high_heterozygosity" --mode=+ - | \
        vk filter ALT --max=0.99 - | \
        vk filter REF --min=1 - | \
        vk filter ALT --min=1 - | \
        vcffixup - | \\
        bcftools view -O z - > filtered.vcf.gz
        bcftools index -f filtered.vcf.gz
    """
}

filtered_vcf.into { filtered_vcf_snpeff; filtered_vcf_to_clean; filtered_vcf_gtcheck; filtered_vcf_phylo }

process annotate_vcf_snpeff {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("merged.filtered.vcf.gz"), file("merged.filtered.vcf.gz.csi") from filtered_vcf_snpeff

    output:
        set file("snpeff.vcf.gz"), file("snpeff.vcf.gz.csi") into snpeff_vcf

    """
        bcftools view -O v merged.filtered.vcf.gz | \\
        snpEff eff ${annotation_reference} | \\
        bcftools view -O z > snpeff.vcf.gz
        bcftools index snpeff.vcf.gz
    """

}


process generate_clean_vcf {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("filtered.vcf.gz"), file("filtered.vcf.gz.csi") from filtered_vcf_to_clean

    output:
        set file("clean.vcf.gz"), file("clean.vcf.gz.csi") into clean_vcf_to_impute
        set val('clean'), file("clean.vcf.gz"), file("clean.vcf.gz.csi") into clean_vcf_stat

    """
        # Generate clean vcf
        bcftools view -m 2 -M 2 --types snps filtered.vcf.gz | \\
        bcftools filter --set-GTs . --exclude 'FORMAT/FT != "PASS"' | \\
        vk filter MISSING --max=0.90 - | \\
        vk filter HET --max=0.10 - | \\
        vk filter REF --min=1 - | \\
        vk filter ALT --min=1 - | \\
        bcftools view -O z > clean.vcf.gz
        bcftools index -f clean.vcf.gz
    """
}


process imputation {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    cpus variant_cores

    input:
        set file("filtered.vcf.gz"), file("filtered.vcf.gz.csi") from clean_vcf_to_impute 
    output:
        set file("impute.vcf.gz"), file("impute.vcf.gz.csi") into impute_vcf
        set val('impute'), file("impute.vcf.gz"), file("impute.vcf.gz.csi") into impute_vcf_stat

    """
        java -jar ${beagle_location} nthreads=${variant_cores} window=8000 overlap=3000 impute=true ne=17500 gt=filtered.vcf.gz out=impute
        bcftools index -f impute.vcf.gz
    """
}

impute_vcf.into { kinship_vcf;  mapping_vcf; }

process make_kinship {

    publishDir analysis_dir + "/cegwas", mode: 'copy'

    input:
        set file("impute.vcf.gz"), file("impute.vcf.gz.csi") from kinship_vcf
    output:
        file("impute.${date}.vcf.gz")

    """
        Rscript -e 'library(cegwas); kinship <- generate_kinship("impute.vcf.gz"); save(kinship, file = "kinship.${date}.Rda");'
    """

}

process make_mapping {

    publishDir analysis_dir + "/cegwas", mode: 'copy'

    input:
        set file("impute.vcf.gz"), file("impute.vcf.gz.csi") from mapping_vcf
    output:
        file("snps.${date}.Rda")

    """
        Rscript -e 'library(cegwas); snps <- generate_mapping("impute.vcf.gz"); save(snps, file = "snps.${date}.Rda");'
    """

}


process calculate_gtcheck {

    publishDir analysis_dir + "/concordance", mode: 'copy'

    input:
        set file("merged.filtered.snp.vcf.gz"), file("merged.filtered.snp.vcf.gz.csi") from filtered_vcf_gtcheck

    output:
        file("gtcheck.tsv") into gtcheck

    """
        echo -e "discordance\\tsites\\tavg_min_depth\\ti\\tj" > gtcheck.tsv
        bcftools gtcheck -H -G 1 merged.filtered.snp.vcf.gz | egrep '^CN' | cut -f 2-6 >> gtcheck.tsv
    """

}

/* 
    Stat VCFs
*/


vcf_stat_set = filtered_vcf_stat.concat( clean_vcf_stat, impute_vcf_stat)

process stat_tsv {

    tag { vcf_stat }

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set val(vcf_stat), file("${vcf_stat}.vcf.gz"), file("${vcf_stat}.vcf.gz.csi") from vcf_stat_set

    output:
        file("${vcf_stat}.stats.txt") into filtered_stats

    """
        bcftools stats --verbose ${vcf_stat}.vcf.gz > ${vcf_stat}.stats.txt
    """

}

/*
    Perform concordance analysis
*/

concordance_script = Channel.fromPath("process_concordance.R")

process process_concordance_results {


    publishDir analysis_dir + "/concordance", mode: "copy"

    input:
        file "gtcheck.tsv" from gtcheck
        file "filtered.stats.txt" from filtered_stats
        file "SM_coverage.tsv" from SM_coverage_merged
        file 'process_concordance.R' from concordance_script

    output:
        file("concordance.svg")
        file("concordance.png")
        file("xconcordance.svg")
        file("xconcordance.png")
        file("isotype_groups.tsv")
        file("gtcheck.tsv") into gtcheck_network

    """
    Rscript --vanilla process_concordance.R
    """

}


filtered_vcf_phylo_contig = filtered_vcf_phylo.spread(["I", "II", "III", "IV", "V", "X", "MtDNA", "genome"])

/*
    Phylo analysis
*/
process phylo_analysis {

    publishDir analysis_dir + "/phylo", mode: "copy"

    tag { contig }

    input:
        set file("merged.filtered.vcf.gz"), file("merged.filtered.vcf.gz.csi"), val(contig) from filtered_vcf_phylo_contig

    output:
        set val(contig), file("${contig}.tree") into trees

    """
        if [ "${contig}" == "genome" ]
        then
            vk phylo tree nj merged.filtered.vcf.gz > genome.tree
        else
            vk phylo tree nj merged.filtered.vcf.gz ${contig} > ${contig}.tree
        fi
    """
}

trees_phylo = trees.spread(Channel.fromPath("process_trees.R"))

process plot_trees {

    publishDir analysis_dir + "/phylo", mode: "copy"

    tag { contig }

    input:
        set val(contig), file("${contig}.tree"), file("process_trees.R") from trees_phylo

    output:
        file("${contig}.svg")
        file("${contig}.png")

    """
    Rscript --vanilla process_trees.R ${contig}
    """

}

process download_annotation_files {
    executor 'local'

    output:
        set val("phastcons"), file("elegans.phastcons.wib") into phastcons
        set val("phylop"), file("elegans.phylop.wib") into phylop
        set val("repeatmasker"), file("elegans_repeatmasker.bb") into repeatmasker

    """
        wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS258/MULTI_SPECIES/hub/elegans/elegans.phastcons.wib
        wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS258/MULTI_SPECIES/hub/elegans/elegans.phylop.wib
        wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS258/MULTI_SPECIES/hub/elegans/elegans_repeatmasker.bb
    """
}

phastcons.mix(phylop).into { wig }

process wig_to_bed {

    publishDir 'tracks', mode: 'copy'

    input:
        set val(track_name), file("track.wib") from wig
    output:
        file("${track_name}.bed.gz") into bed_tracks
        file("${track_name}.bed.gz.tbi") into bed_indices

    """
        bigWigToBedGraph track.wib ${track_name}.bed
        bgzip ${track_name}.bed
        tabix ${track_name}.bed.gz
    """

}

process final_vcf {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    cpus alignment_cores

    input:
        set file('snpeff.vcf.gz'), file('snpeff.vcf.gz.csi') from snpeff_vcf
        file(track) from bed_tracks.toSortedList()
        file(track) from bed_indices.toSortedList()
        file('vcf_anno.conf') from Channel.fromPath("vcfanno.conf")

    output:
        set file("WI.${date}.vcf.gz"), file("WI.${date}.vcf.gz.csi"), file("WI.${date}.vcf.gz.tbi")
    """
        vcfanno -p ${alignment_cores} vcf_anno.conf snpeff.vcf.gz | bcftools view -O z > WI.${date}.vcf.gz
        bcftools index WI.${date}.vcf.gz
        tabix WI.${date}.vcf.gz
    """

}


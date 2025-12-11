#####bash script#####

####vcftools filtering
#snps only

vcftools --vcf populations.snps.vcf --remove lowDP.indv --max-missing 0.75 --max-alleles 2 --minDP 3 --maf 0.05 --min-meanDP 20 --recode --recode-INFO-all --out feems_filt

perl -pe 's/^([^#])/scaffold_\1/' feems_filt.recode.vcf > feems_filt.recode.chr.vcf

plink --vcf feems_filt.recode.chr.vcf --make-bed --out feems_bed --allow-extra-chr

VCF2PCACluster -InVCF feems_filt.recode.vcf -OutPut pca_feems -InSampleGroup ~/shovel-bugs/new_pop_file.txt 

#invariant + variant file for pixy

vcftools --vcf populations.all.vcf --remove lowDP.indv --recode --recode-INFO-all --out populations_no_lowDp
vcftools --vcf populations_no_lowDp.recode.vcf --remove ~/shovel-bugs/cryptic_species.txt --recode --stdout | bgzip -c > initial_filtered_indv.vcf.gz 

#separate the invariants from variants

vcftools --gzvcf initial_filtered_indv.vcf.gz --remove-indels --max-alleles 1 --recode --stdout | bgzip -c > invariant_sites.vcf.gz #this keeps invariants

vcftools --gzvcf invariant_sites.vcf.gz --max-missing 0.75 --minDP 20 --recode --stdout | bgzip -c > test_missing_invariant_sites.vcf.gz #this keeps fewer invariants

# create a filtered VCF containing only variant sites
vcftools --gzvcf initial_filtered_indv.vcf.gz --mac 1 --remove-indels --max-missing 0.75 --maf 0.05 --min-meanDP 20  --recode --stdout | bgzip -c > variant_sites.vcf.gz

# index both vcfs using tabix
tabix test_missing_invariant_sites.vcf.gz
tabix variant_sites.vcf.gz

# combine the two VCFs using bcftools concat
bcftools concat --allow-overlaps test_missing_invariant_sites.vcf.gz variant_sites.vcf.gz -O z -o final_filtered_sites.vcf.gz

bcftools view --header-only final_filtered_sites.vcf.gz | grep "##contig" | cut -f3 -d "=" | sed 's/>/''/g' | awk '{print $0, " chr" int((NR-1)/10000+1)}'  >chrlist.txt

bcftools annotate --rename-chrs chrlist.txt final_filtered_sites.vcf.gz -Oz -o final_renamed.vcf.gz

bcftools view -H final_renamed.vcf.gz > rename.body.vcf

zgrep -v "^#" final_renamed.vcf.gz |wc -l
echo {1..10603559} | tr ' ' '\n' > contpos.txt

awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[NR]=$1; next} {if (FNR in a) $2=a[FNR]; print $0}' contpos.txt rename.body.vcf > rename2.body.vcf

bcftools sort final_renamed.vcf.gz -Oz -o final_sorted.vcf.gz

bcftools view -h final_sorted.vcf.gz > sort.head.txt

cat sort.head.txt rename2.body.vcf | bgzip > final_merged.vcf.gz

tabix final_merged.vcf.gz

##pixy analysis

pixy --stats fst \
--vcf final_merged.vcf.gz \
--populations new_pop_file.txt \
--window_size 10603559 \
--n_cores 4

pixy --stats pi \
--vcf final_merged.vcf.gz \
--populations new_indv_file.txt \
--window_size 10603559 \
--n_cores 4

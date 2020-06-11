######## vcf to bfile ######
vcf="/home/lhjiang/data/arabidopsis/1001genome/1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf"
plink --vcf $vcf --make-bed --out arabi # vcf to bed
cp arabi.bim arabi.bim.backup
awk '{print $1, $1"_"$4, $3, $4, $5, $6}' arabi.bim.backup > arabi.bim
plink --bfile arabi --recode vcf-iid --out arabi
######### vcf impute ########
# left
awk '(NR>11){OFS="\t";print $1, $2, $3, $4, $5, $6, $7, $8, $9}' arabi.vcf > arabi.left 
# right test
#awk '{a=""; for(i=1;i<=NF;i++){a=a$i"/"$i"\t"}; print a}' data/snps.10.txt > data/arabi.10.right
# right
awk '{a=""; for(i=1;i<=NF;i++){a=a$i"/"$i"\t"}; print a}' snps.txt > arabi.right
# head
head -n 11 arabi.vcf > arabi.head
# merge left and right 
paste arabi.left arabi.right > arabi.imputed.nohead.vcf
# add head
cat arabi.head  arabi.imputed.nohead.vcf > arabi.imputed.vcf
plink --vcf arabi.imputed.vcf --make-bed --out arabi.imputed

#awk '{print  $1"_"$4}' atwell.bim > atwell.snps
plink --bfile arabi.imputed --geno 0.05 --extract atwell.snps --make-bed --keep keep.txt --pheno FT10.plink --out arabi

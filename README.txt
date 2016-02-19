
* TRAP 2.3 - README *

###### The structure of directories ######
  ./tophat_results/ 	: where tophat results are stored
	/SAMPLE1/
	/SAMPLE2/
	/ ...  /
	/SAMPLEN/
  ./cufflinks_results/	: where cufflinks results are stored
	/SAMPLE1/
	/SAMPLE2/
        / ...  /
        /SAMPLEN/
  ./cuffdiff_results/ 	: where cuffdiff results are stored
	/DIFF1/
	/ ... /
	/DIFFN/
  ./oryza_sativa/	: where the files used for rice analysis are stored (gff, kgml, ...)
  ./TRAP_results/            : where analysis results will be stored


###### How to run the program ######
1. If you already have Tophat, Cuffdiff and/or Cufflinks result, put the output files into each sample directory.
   (If you have Cufflinks or Cuffdiff result, TRAP does not need Tophat results.)
  * Essential file for Tophat :		accepted_hits.bam
  * Essential file for Cufflinks : 	genes.fpkm_tracking
  * Essential file for Cuffdiff : 	gene_exp.diff	

2. Open and modify 'config.txt' with your own configuration.
   There are some optional input which is nesessary for Tophat and Cufflinks.
 
3. Run SPIA (type ./TRAP)


###### About configuration in config.txt ######
As for 3, 4, 5, We provide the files as default for Oryza Sativa L. ssp. japonica)
All the path should be written from the root. 
(Eg. TRAP/oryza_sativa (X), /PATH_TO_TRAP/TRAP/oryza_sativa (O))

Assume you have RNA-seq data from 2 normal and 2 cancer patients, 
sequenced at three different time points (10 day, 20 days, 30 days after applying drug).
1. Number of time points
  numTP=3
2. Sample pair : User-defined sample names. 
  If you are to use sample name CON for control, TRE for treatment, write :
  control1=CON10
  control2=CON20
  control3=CON30
  treatment1=TRE10
  treatment2=TRE20
  treatment3=TRE30

  If there are more than two replicates for each sample, write each directory
name separated by comma as :
  control1=CON11,CON12,CON13

3. Gene name conversion file path
  This is for conversion between KEGG gene IDs and gene symbols in Cufflinks result.
  Gene symbols are read from 'gene_short_name' column in the 'genes.fpkm_tracking'
file.

  1) Make your own tab-delimited file as follows :
        hsa:10186       LHFP
        hsa:10814       CPLX2
        ...

  2) Or download file from http://rest.kegg.jp/list/***
     (*** is an organism code from http://www.kegg.jp/kegg/catalog/org_list.html)
     If multiple gene symbols are matched with one KEGG gene ID, TRAP uses the first one.
     >>>>>>>>> Use this when the gene symbols in the file is same as those in Cufflinks result.

4. Pathay name file path
  Download from http://rest.kegg.jp/list/pathway/***	
  (*** is an organism code from http://www.kegg.jp/kegg/catalog/org_list.html)

5. KGML file path
  Path to the Directory containing kgml files downloaded from http://rest.kegg.jp/get/***00000/kgml
  (***00000 is a pathway code from http://rest.kegg.jp/list/pathway/***)
  File names should be '***00000.xml'.
  
6. Use of Cuffdiff for finding DEGs
  If you want to use Cuffdiff results for finding DEGs by p-value,
  write 'yes' and specify the cutoff value for p-value.

6-1. Cuffdiff directory name
  Write only if you answered 'yes' in 6.
  Write down the directory names for each time point.
  Ex. Assume you compared CON10 and TRE10, 
      and the results are in cuffdiff_results/DIFF10/gene_exp.diff. 
      Then write :
  diff1=DIFF10
  diff2=DIFF20
  ...

6-2. DEG/cluster cutoff (log fold change)
  Write only if you answered 'no' in 6.

7. Time-lag factor
  The ratio of the interaction from the upstream genes of current vs. previous time point.
  If downstream genes are affected only by the previous time point, set as 1.0.

8. Tophat software path
9. Cufflinks software path
  Path to Tophat/Cufflinks in your server.

10. Sample name and file path 	(optional for Tophat)
  Left side of the equal sign should be "sample name you defined" + "Path"
  Paired-end sequencing data files should be seperated by comma.
  (Eg. CON10Path=/PATH_TO_THE_FILE/control10_1.fq,/PATH_TO_THE_FILE/control10_2.fq)

11. Reference genome index path (optional for Tophat) 
  If you indexed reference genome by Bowtie, several files (Eg. REFERENCE_FILE_NAME.X.bt) are generated.
  Write the path to the files including reference genome file name.
  (Eg. refIndex=/PATH_TO_THE_FILE/REFERENCE_FILE_NAME )

12. GFF file path 		(optional for Cufflinks)
  Path to the supplied reference annotation.
  This is for -G option in Cufflinks. 
  The gene names in the annotation file should be same as those in the gene conversion file (3).
  

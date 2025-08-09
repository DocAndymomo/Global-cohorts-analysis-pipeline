{\rtf1\ansi\ansicpg936\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 Times-Bold;\f1\froman\fcharset0 Times-Roman;\f2\froman\fcharset0 Times-Italic;
\f3\fmodern\fcharset0 Courier-Bold;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red109\green109\blue109;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\cssrgb\c50196\c50196\c50196;}
{\*\listtable{\list\listtemplateid1\listhybrid{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{disc\}}{\leveltext\leveltemplateid1\'01\uc0\u8226 ;}{\levelnumbers;}\fi-360\li720\lin720 }{\listname ;}\listid1}
{\list\listtemplateid2\listhybrid{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{disc\}}{\leveltext\leveltemplateid101\'01\uc0\u8226 ;}{\levelnumbers;}\fi-360\li720\lin720 }{\listname ;}\listid2}}
{\*\listoverridetable{\listoverride\listid1\listoverridecount0\ls1}{\listoverride\listid2\listoverridecount0\ls2}}
\paperw11900\paperh16840\margl1440\margr1440\vieww18760\viewh21340\viewkind0
\deftab720
\pard\pardeftab720\sa321\partightenfactor0

\f0\b\fs48 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Global Cohorts Analysis Pipeline\
\pard\pardeftab720\sa298\partightenfactor0

\fs36 \cf0 Introduction\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 The 
\f0\b Global Cohorts Analysis Pipeline
\f1\b0  enables the 
\f0\b systematic recovery
\f1\b0  and 
\f0\b high-precision profiling
\f1\b0  of microbial communities from human genomic and transcriptomic sequencing data, with 
\f0\b minimized contamination
\f1\b0 .\
This pipeline integrates 
\f0\b two-step host sequence removal
\f1\b0 , 
\f0\b robust species annotation
\f1\b0 , and 
\f0\b multi-center contamination filtering
\f1\b0 , making it highly suitable for 
\f0\b large-scale, cross-cohort microbiome studies
\f1\b0 .\
By combining 
\f0\b deep host depletion
\f1\b0 , 
\f0\b high-resolution taxonomic classification
\f1\b0 , and 
\f0\b multi-cohort contamination filtering
\f1\b0 , the pipeline achieves 
\f0\b lower contamination rates
\f1\b0  and 
\f0\b greater analytical precision
\f1\b0  than conventional metagenomic workflows.\cf3 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Overview\
\pard\pardeftab720\sa240\partightenfactor0

\fs24 \cf0 Main Features
\f1\b0 \
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls1\ilvl0
\f0\b \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Two-step host removal
\f1\b0  using both 
\f2\i GRCh38.p14
\f1\i0  and 
\f2\i CHM13 v2.0
\f1\i0  reference genomes to maximize removal of host sequences.\
\ls1\ilvl0
\f0\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 High-accuracy species annotation
\f1\b0  with 
\f2\i Kraken2
\f1\i0  + 
\f2\i Bracken
\f1\i0  using the 
\f2\i PlusPF
\f1\i0  database.\
\ls1\ilvl0
\f0\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Cross-cohort contamination removal
\f1\b0  using a 
\f0\b multi-center contamination database
\f1\b0  and 
\f2\i decontam
\f1\i0  filtering.\
\ls1\ilvl0
\f0\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Final microbiome abundance tables
\f1\b0  generated from high-confidence reads for downstream statistical analyses.\
\pard\pardeftab720\partightenfactor0
\cf3 \strokec3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Step-by-Step Pipeline\
\pard\pardeftab720\sa280\partightenfactor0

\fs28 \cf0 Step 1 \'96 Raw Data Conversion (
\f3\fs30\fsmilli15210 01-Raw2fastq
\f0\fs28 )\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 Convert raw sequencing output (e.g., 
\f0\b BCL format
\f1\b0 ) into 
\f0\b FASTQ
\f1\b0  files, ensuring proper sample demultiplexing and metadata integrity.\
\pard\pardeftab720\sa280\partightenfactor0

\f0\b\fs28 \cf0 Step 2 \'96 Quality Control (
\f3\fs30\fsmilli15210 02-fastp
\f0\fs28 )\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 Perform read trimming and quality filtering using 
\f0\b fastp (v0.23.4)
\f1\b0  to remove low-quality bases, adapters, and sequencing artifacts.\
\pard\pardeftab720\sa280\partightenfactor0

\f0\b\fs28 \cf0 Step 3 \'96 First Host Removal (GRCh38) (
\f3\fs30\fsmilli15210 03-bowtie2_grch
\f0\fs28 )\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 Align reads to the 
\f0\b human reference genome GRCh38.p14
\f1\b0  using 
\f0\b bowtie2 (v2.5.2)
\f1\b0  to remove human-origin sequences.\
\pard\pardeftab720\sa280\partightenfactor0

\f0\b\fs28 \cf0 Step 4 \'96 Second Host Removal (CHM13) (
\f3\fs30\fsmilli15210 04-bowtie2_chm
\f0\fs28 )\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 Align the remaining reads to the 
\f0\b CHM13 v2.0
\f1\b0  telomere-to-telomere (T2T) human reference using 
\f0\b bowtie2 (v2.5.2)
\f1\b0  for deeper host sequence depletion.\
\pard\pardeftab720\sa280\partightenfactor0

\f0\b\fs28 \cf0 Step 5 \'96 Taxonomic Classification & Quantitation (
\f3\fs30\fsmilli15210 05-kraken_beacken
\f0\fs28 )\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 Classify microbial reads using 
\f0\b Kraken2 (v2.1.3)
\f1\b0  with the 
\f0\b PlusPF (2024.1.12)
\f1\b0  database.\uc0\u8232 Refine species-level abundance estimates using 
\f0\b Bracken
\f1\b0  for more accurate microbial quantitation.\
\pard\pardeftab720\sa280\partightenfactor0

\f0\b\fs28 \cf0 Step 6 \'96 Contamination Removal (
\f3\fs30\fsmilli15210 06-Decontam
\f0\fs28 )\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 Use a 
\f0\b multi-center contamination database
\f1\b0  and 
\f0\b decontam
\f1\b0  filtering to remove background contaminants, ensuring cross-cohort comparability.\
\pard\pardeftab720\sa280\partightenfactor0

\f0\b\fs28 \cf0 Step 7 \'96 Microbiome Data Extraction (
\f3\fs30\fsmilli15210 07-kraken_tools
\f0\fs28 )\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 Extract high-confidence microbial reads and generate final 
\f0\b microbiome abundance tables
\f1\b0  using 
\f0\b KrakenTools
\f1\b0 , ready for downstream statistical and ecological analyses.\
\pard\pardeftab720\partightenfactor0
\cf3 \strokec3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Advantages Over Conventional Methods\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls2\ilvl0
\fs24 \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Reduced Contamination
\f1\b0  \'97 Multi-step host removal and cohort-specific contaminant filtering.\
\ls2\ilvl0
\f0\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Improved Precision
\f1\b0  \'97 Dual-host genome alignment combined with Bracken refinement for species-level accuracy.\
\ls2\ilvl0
\f0\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Scalable
\f1\b0  \'97 Optimized for integration of datasets from multiple global cohorts.\
\ls2\ilvl0
\f0\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Reproducible
\f1\b0  \'97 Modular structure with version-controlled tools and databases.\
}
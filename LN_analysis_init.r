# This file has includes and static data used across LN analysis
# Run this script prior to any LN_analysis*.r files, cell_cycle_regression.r,
# and volcano_plot.r.

library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)

cytotoxic_gene_list <- c("Prf1", "Gzmb", "Gzmk", "Ccl4", "Ccl5", "Csf1")
traffic_gene_list <- c("S1pr1", "Sell", "Ccr7", "Cxcr4", "Cxcr3", "Cxcr6")
gen_gene_list <- c("Itgae", "Ccr7", "Klf2", "Cxcr6", "S1pr1", 
                  "Cd8a", "Thy1", "Ptprc", "Cd3e")
custom_color_palette <-  c("#1B9E77", "#D95F02", "#7570B3", "#CA00C4",
                          "#66A61E", "#E6AB02", "#A6761D", "#666666")


LN_Trm_genes <- c("Acap1", "Actn2", "Amica1", "Arhgef1", "Atxn7l3b", "Aw112010",
                  "B4galnt1", "Bcl11b", "Cbx3", "Ccnd2", "Ccr10", "Cd27", "Cd7",
                  "Cd74", "Chd3", "Cirbp", "Clec2d", "Crot", "Csf1", "Cxcr3",
                  "Cxcr6", "Eif5", "Evl", "Fam189b", "Fubp1", "Fyb", "Gramd1a",
                  "Sema4a","Gstp1", "Shisa5", "H2-T23", "Sipa1", "Hmgn1",
                  "Slfn2", "Hmha1", "Sp100", "Hsp90b1", "Spcs2", "Id2", "Srrm2",
                  "Ifitm10", "Stap1", "Ikzf3", "Tbc1d10c", "Il16", "Tesc",
                  "Il18r1", "Tnfaip8", "Il7r", "Tnrc6a", "Irf2bpl", "Tsc22d4",
                  "Itgae", "Uba52", "Itgal", "Ucp2", "Itm2c", "Wbp1", "Lfng",
                  "Xist", "Lpar6", "Ypel3", "Lrrc58", "Zbtb7a", "Ltb", "Znrf1",
                  "Ly6a", "Ly6e", "Ly6g5b", "Malat1", "Mbnl1", "Mrpl52", "Mxd4",
                  "Mycbp2", "N4bp2l2", "N4bp2l2", "Ndufa3", "Ndufa5", "Nktr",
                  "Nudcd3", "Ogt", "Pdcd4", "Pdia3", "Pdia6", "Ptpn7", "Ptprc",
                  "Rapgef6", "Rbpj", "Rgs10", "Rpl15", "Rpl35", "Rpl38",
                  "Rps28", "Rps29", "Sash3")
Tcm_genes <- c("1500002C15Rik","1700029J07Rik","1700041G16Rik","1700058P15Rik",
               "1700061G19Rik","1700096K18Rik","2610027K06Rik","2810405F17Rik",
               "3110082I17Rik","4930417O13Rik","4930522L14Rik","4930579G18Rik",
               "5830462O15Rik","Abcd3","AC100406.1","AC117769.3","AC124556.1",
               "AC125409.1","AC126549.2","AC153959.1","AC154377.1","AC157515.1",
               "AC160637.1","AC167565.1","Actn1","Aff3","AI606181","Ankrd37",
               "Anxa11os","Anxa9","Apbb2","Apoe","Aqp11","Arc","Arhgap21",
               "Arhgef10","Arhgef39","Art4","Atp6v1g3","B430306N03Rik",
               "BC055324","Bend4","Bhlha15","Birc5","Bivm","Blzf1","Bmf",
               "Cacnb1","Car12","Ccdc142os","Ccdc162","Ccdc78","Ccnb2",
               "Ccnb2-ps","Ccr7","Cdca2","Cdca7","Cdr2l","Celsr2","Cenpk",
               "Cenpu","Cep41","Ckap2","Clcn2","Clcn5","Cldn12","Cldnd2",
               "Cnnm2","Copg2","Cpq","Cstf2t","CT025639.2","Ctxn1","Cxcr5",
               "Cyb5d2","Cym","Cyp2u1","D030055H07Rik","D130020L05Rik",
               "D730003I15Rik","D830050J10Rik","Dapl1","Dctd","Degs2","Dglucy",
               "Dhdh","Dnah11","Doc2g","Dsn1","E030042O20Rik","E130215H24Rik",
               "E2f3","E330034L11Rik","Egf","Eomes","Erf","Erp27","Ets2",
               "Fam57a","Fam71e1","Fmo5","Frg2f1","Galc","Galnt4","Gcsh",
               "Gfod2","Gipc3","Glb1l","Glce","Gm11263","Gm11342","Gm12258",
               "Gm12596","Gm12743","Gm12912","Gm13708","Gm15441","Gm15470",
               "Gm15503","Gm15550","Gm15559","Gm15672","Gm15788","Gm15903",
               "Gm15953","Gm16150","Gm16152","Gm16209","Gm17195","Gm17259",
               "Gm19589","Gm20707","Gm26648","Gm26852","Gm36975","Gm37065",
               "Gm37069","Gm37137","Gm37305","Gm37422","Gm37537","Gm37570",
               "Gm38034","Gm38102","Gm38248","Gm38336","Gm42432","Gm42572",
               "Gm43088","Gm43167","Gm43292","Gm43339","Gm43340","Gm43581",
               "Gm43654","Gm43737","Gm44510","Gm44616","Gm44888","Gm45051",
               "Gm45290","Gm45406","Gm45902","Gm7434","Gm8325","Gm9645",
               "Gm9917","Gpr15","Gtdc1","Has3","Hcn3","Hells","Hist1h2be",
               "Hk2","Ifngas1","Igsf3","Il6ra","Iqcb1","Katnal1","Kcnip3",
               "Kcnk7","Kifap3","Kifc3","Klhl7","Lancl2","Lbhd1","Lclat1",
               "Lipt2","Lrig1","Lrrc23","Mapk11","Mapk8ip2","Mapkapk5","Mast1",
               "Mast2","Mboat1","Mboat4","Mcm8","Mid1","Mki67","Mmp11","Mpp5",
               "Mpv17l","Mr1","Mtbp","Mtfp1","Nacc2","Ncf2","Ndrg2","Neu3",
               "Neurl2","Nfxl1","Nipal1","Nle1","Nr4a1","Nsg2","Nt5e","Oasl2",
               "Paqr8","Pdzd2","Pfkfb1","Pfkfb2","Piga","Pigg","Pkd1l3","Pld4",
               "Pnpla8","Poc1a","Polr1b","Pomgnt1","Ppic","Prkag2","Prss2",
               "Prss36","Pygm","Rad51","Rapgef4","Rasgef1b","Rasgef1c","Rcan1",
               "Rec114","Rft1","Rgs11","Rph3al","Rpp40","Rps19-ps5","Rps4l",
               "Rwdd3","Scarf1","Sdhaf3","Sell","Sfxn4","Shmt1","Siah2",
               "Sipa1l2","Skp2","Slc12a4","Slc16a5","Slc25a13","Slc29a2",
               "Slc37a2","Slc43a1","Slc6a19","Smim13","Smyd2","Socs6","Specc1",
               "Spef1","Spef2","Ss18l1","Stxbp4","Susd1","Tbc1d8b","Tceanc",
               "Tctex1d4","Tctn1","Tdrp","Thnsl1","Thnsl2","Tjp3","Tlr1",
               "Tmem14a","Tmem170","Tmem231","Tmem50b","Tmem51","Tmie",
               "Tnfaip8l1","Tnfrsf12a","Tnfrsf13c","Tns1","Top1mt","Tppp3",
               "Treml2","Trim32","Trmt10a","Trp53inp1","Ttc16","Tti1","Tvp23a",
               "Tymp","Ube2c","Usp2","Vpreb1","Vwa7","Wdr25","Wdr62","Wnk4",
               "Xk","Zfp189","Zfp202","Zfp296","Zfp41","Zfp438","Zfp516",
               "Zfp568","Zfp583","Zfp773","Zkscan5")
Tem_genes <- c("0610043K17Rik","1500009L16Rik","1700120C14Rik","1810011H11Rik",
               "1810062G17Rik","2010109I03Rik","2810468N07Rik","4930509H03Rik",
               "4930599N23Rik","6330403K07Rik","9030407P20Rik","9130221H12Rik",
               "9130604C24Rik","9230112E08Rik","9330198N18Rik","A430035B10Rik",
               "A530083M17Rik","A930037H05Rik","AA467197","Abca3","Abhd13",
               "AC093043.2","AC115863.1","AC116680.1","AC122197.1","AC122273.1",
               "AC153506.1","AC163354.1","Acbd6","Acer2","Acot11","Ada",
               "Adam15","Adam8","Adamts14","Adap1","Adssl1","Afap1","Ager",
               "Agrn","AI846148","Akap12","Alcam","Aldh7a1","Alyref","Ankrd29",
               "Ano8","Anxa4","Aplf","Aplp1","Apobec2","Apol7b","Arid3a",
               "Arl10","Art2a-ps","Asb2","Asb7","Asns","Asph","Aspm","Atp2a1",
               "Atp6v1g2","Atp8a2","Atxn7l1os2","B230219D22Rik","B3gat2",
               "BC052040","BC064078","Bcl2a1d","Bclaf3","Bend6","Bhlhe40","Bik",
               "Bloc1s3","Bloc1s4","C030034I22Rik","C030034L19Rik","Cacna1c",
               "Cage1","Calr3","Car7","Caskin2","Ccdc102a","Ccdc136","Ccdc141",
               "Ccr9","Ccrl2","Cd160","Cd200r1","Cd22","Cd86","Cd93","Cdkl2",
               "Cdkn1a","Cdkn2c","Cep170b","Cfap126","Chn2","Chst14","Chst15",
               "Cish","Cited2","Ckb","Col9a3","Coro2a","Cox11","Crb2","Crybg2",
               "Csgalnact2","Csnk1e","Csrnp2","CT025642.2","Ctnnbip1","Ctsf",
               "Cx3cr1","Cxcl10","Cyp2s1","Cyp4f16","D030056L22Rik",
               "D3Ertd751e","D630030B08Rik","Dapk3","Ddx28","Dgat2","Dhrs13",
               "Dnajb5","Dnajb7","Dnajc12","Dnajc28","Dock6","Dusp2","E2f2",
               "E2f5","Efemp2","Eid2b","Eme1","Eng","Enkd1","Enpp5","Epcam",
               "Ephx1","Fads1","Fah","Fam122a","Fam129b","Fam20a","Fam212a",
               "Fancf","Fbxl15","Fer1l5","Fes","Fgl2","Fhl3","Flywch2","Frmd4b",
               "Galnt3","Gch1","Gclm","Gcnt1","Gdpd1","Gem","Ggt5","Gm1043",
               "Gm10522","Gm10719","Gm10863","Gm11702","Gm12352","Gm12542",
               "Gm13349","Gm13623","Gm14125","Gm15333","Gm15337","Gm15518",
               "Gm15601","Gm15787","Gm16364","Gm16576","Gm16702","Gm17477",
               "Gm17764","Gm2011","Gm20257","Gm20511","Gm20628","Gm26614",
               "Gm26799","Gm2788","Gm28100","Gm28530","Gm31597","Gm31763",
               "Gm34220","Gm36981","Gm37663","Gm37677","Gm38091","Gm38227",
               "Gm39041","Gm4117","Gm4208","Gm42372","Gm42571","Gm42576",
               "Gm42735","Gm42819","Gm43361","Gm43637","Gm43857","Gm43963",
               "Gm44292","Gm4430","Gm44981","Gm45548","Gm45592","Gm45698",
               "Gm45779","Gm6637","Gm9752","Gnaq","Gpam","Gpatch3","Gpnmb",
               "Gpr34","Gpr55","Grasp","Grb7","Grin3b","Gstt1","Gzma","Gzmb",
               "H1f0","Hic1","Hist2h2ac","Hmox1","Hspa2","I830127L07Rik",
               "Icosl","Ift88","Il1rl1","Il2ra","Iqcg","Irf5","Itga1","Itih5",
               "Itpkc","Itsn1","Jmjd7","Kcnj8","Kctd12","Kctd17","Kctd7",
               "Kdm1b","Khdrbs1","Kif5c","Klc3","Klhdc8b","Klrb1c","Klre1",
               "Klrg1","Krt83","Lactbl1","Lag3","Lamc1","Lamp1","Lgals3",
               "Lhpp","Lilr4b","Lilrb4a","Litaf","Lpar2","Lrp5","Lrrc27",
               "Lrrc4","Lzts1","Macrod1","Maneal","Mansc1","Map3k15","Matk",
               "Med16","Mfsd13a","Mgat4b","Miat","Mical2","Mpi","Mpp6","Msh5",
               "Msh6","Mt1","Muc1","Mxi1","Myo6","Napsa","Ndc1","Nebl","Nek6",
               "Nhlrc1","Nhsl2","Nmb","Nmral1","Noct","Nrg4","Nsmce3","Ogfrl1",
               "Optn","Orai1","Osbpl3","P2rx7","Palm","Pcbp4","Pcgf6",
               "Pdcd1lg2","Pdcd7","Pde11a","Pde4a","Pdlim7","Perm1","Pfkfb4",
               "Pfn2","Pgam2","Phldb3","Pisd-ps2","Plcxd1","Plekhf1","Plk1",
               "Plxna1","Plxna3","Plxnb1","Ppfibp1","Ppm1d","Ppp1r2","Ppp1r35",
               "Ppp1r3e","Prdm16","Prkra","Prrt1","Ptgdr2","Ptges2","Ptms",
               "Ptov1","Qpct","Rab11fip2","Rab24","Rab43","Rac1","Ralb",
               "Rasl11a","Rbbp9","Rd3","Rdh10","Rdh5","Rgs1","Rgs12","Rhoq",
               "Ripor3","Rmdn2","Rmi2","Rnase4","Rnf187","Rnf216","Rnf39","Rp2",
               "Rpl13-ps5","Rpl23a-ps2","Ryr1","S1pr5","Sap30","Sart3","Sbf2",
               "Scand1","Selenon","Selenoo","8-Sep","Serpina3f","Serpina3g",
               "Serpinc1","Sft2d3","Sh3d19","Slamf1","Slc16a10","Slc18a1",
               "Slc24a1","Slc35e4","Slc39a8","Slc5a11","Slco4a1","Smoc1",
               "Snhg14","Snx20","Socs2","Spata5l1","Spats2","St14","St6galnac3",
               "Stk32c","Stkld1","Sult2b1","Sumf1","Susd2","Syap1","Tbkbp1",
               "Tbxas1","Tcf7l1","Tcp11","Tdrd7","Tef","Tesc","Thap6","Tigd5",
               "Tigit","Timm29","Tmem151a","Tmem160","Tmem176b","Tmem198b",
               "Tmprss13","Traj28","Traj36","Traj60","Traj8","Trat1","Trav8-2",
               "Trim27","Trpt1","Ttyh2","Ube2t","Uhrf1","Vsig10l","Wdr46-ps",
               "Wfikkn1","Wfs1","Wnt2b","Xcr1","Yars2","Zbtb11os1","Zbtb32",
               "Zbtb7b","Zeb2","Zfand5","Zfp28","Zfp341","Zfp358","Zfp467",
               "Zfp551","Zfp580","Zfp683")
TGFb_genes <- c("Mcpt2","Cdh1","Spsb1","Vipr2","Gpr56","Src","Ppp2r2c","Lrig1",
                "Itgae","Agap1","Ncmap","Pmepa1","Sema6d","Emid1","Cd33","Dlk1",
                "Ldlrad4","Car2","Cpd","Nt5e","Tspan9","Gsg2","Klhl30",
                "1810011H11Rik","Osgin1","Ccl1","Litaf","Itga1","Kifc3",
                "Hsf2bp","Asic3","Abi3","Smurf2","Phactr2","Oplah","Qpct",
                "Tfr2","Isg20","Rnase6","Rgs1","2900026A02Rik","Mmp11",
                "Tnfsf11","Nrarp","Cyb561","Smyd1","Kcnip2","Cx3cr1","Nek6",
                "Nlrp1b","St8sia1","Arhgap39","Jup","Htra3","Rgs16","H2-M5",
                "Chn2","Cish","Atp6v0a1","Skil","Dok3","Igflr1","Ccr8","Timp2",
                "Zfyve28","Ppm1n","Hpgds","B4galnt4","Ifng","Ctnnal1","Clec12a",
                "Exoc3l","Coro2a","Ikzf4","Adamts6","D8Ertd82e","Smpd5","Aqp3",
                "Evpl","Ramp1","St8sia6","Xcl1","Scn1b","Rnf149","Dtx4","Gngt2",
                "Sbk1","Tbc1d16","Tnfrsf13c","Gna12","Ermn","Neu3","Fmnl3",
                "Cd83","Epb4.1l2","Ccdc112","Adam19","Rab26","Fam101b","Mical3",
                "Prkcz","Grina","Slc27a6","Tgfbr3","Fgfr1","Msc","Rgs10",
                "Lonrf1","Lax1","Kcnc1","Nphp1","Slc16a10","Kif13a","Ninj1",
                "Smyd3","9430020K01Rik","Csgalnact1","Gpaa1","Ski","Gcnt4",
                "Map9","Egr3","Fam161a","Egr1","Fndc3a","Mapkapk3","Ctss",
                "Hnrnpll","Galm","Dusp2","Stom","Esm1","1700049G17Rik",
                "Plekho1","Med10","Smtn","Gpr34","Sepn1","Egr2","Prrt2","Aen",
                "Cd101","Gtf2ird1","Tiam1","Camkk1","D430042O09Rik","Fam214b",
                "Matk","Ralgps1","Dapk2","Usp6nl","Foxred2","Wdyhv1","Znrf1",
                "Tjp1","Irf8","Hemk1","Pgap1","Accs","Aim2","Per3","Zfr2",
                "Lgalsl","1700001L05Rik","Zfp820","D3Ertd254e","Gcnt1",
                "Slc41a2","Ttc39b","Gclm","Peg13","Slc9a1","Adora3","Cers6",
                "Ccrn4l","Cd96","Golim4","Lpcat2","Lsr","Acsbg1","Eef2k",
                "Plekhf1","Rbm20","Ssx2ip","Ankrd50","Igfbp4","Inpp4b",
                "Irf2bpl","Pygl","Zfp1","Golm1","Gpr68","Ptgfrn","Tsc22d1",
                "Abca1","Fam124b","Itpripl2","Bcl6","Lysmd2","Trp53inp2",
                "Zdhhc13","Bpgm","F2r","Frmd4b","Ctsw","Swap70","Frmd6","Gas7",
                "Gdpd5","Spire1","Tet3","Batf3","Dstyk","Luzp1","Mgat5","Ptpre",
                "Ralgps2","Mif4gd","Stat1","Ttc3","Abhd15","Cerk","Adssl1",
                "Pcyt1b","Rai1","Blcap","Map3k14","Rnf19b","Scai","Tmem57",
                "Atp6v1g2","Chst12","Fam20a","Gtf3c1","Trp53inp1","Wdr78",
                "Aars2","Cd244","Ly6g5b","Tbx6","Usp22","Zfp827",
                "1600014C10Rik","Als2","Arhgef5","B4galt5","Nfat5","Prkacb",
                "Rgs2","Slc9a3r1","Soat1","Tctn3","Ttc39c","Cotl1","Ldlrap1",
                "Ncf1","Iigp1","Ikzf3","Ipcef1","Irf4","Abi2","Runx3","Ypel3",
                "Entpd1","Fut8","Inpp5f","Apol7e","Arhgef12","Nrp1","Slc26a11",
                "Tnfrsf1b","Cd160","Gfod1","Gm12185","H6pd","Pmm1","Tmem2",
                "Ublcp1","Dennd3","Gramd1a","Idh2","Ppip5k1","Slc39a13",
                "Baiap3","Extl3","Mxd1","Nipal1","Rrp1b","Twsg1","Cdc42bpg",
                "Celsr1","Ehd1","Kit","Slc22a15","Tmcc1","Camsap2","Klhl25",
                "Ncf4","Plcxd2","Rab11fip4","Specc1","Fam3c","Fuca2","Pde4a",
                "Prr12","Ctnnb1","Egln3","Fam46a","Fbxo25","Gprin3","Scly",
                "AW112010","Cd1d1","Lrrc61","Clstn1","Exosc4","Smad7","Susd3",
                "Traf4","Vasp","Gne","Gpbp1l1","Prkch","Rab37","Rbpj","Usp11",
                "AA467197","Bmpr2","Cd8a","Dpp9","Inpp5d","Kif1b","Lasp1",
                "Rftn1","Wee1","Fasl","Nbas","Plscr1","Prkdc","Rhoh","Spcs2",
                "Suox","Tbc1d4","Tgif1","Anp32a","Lnpep","Myo5a","Rreb1",
                "Zfp706","Ermp1","Fam149b","Glrx","Pacsin2","Plekha2","Sorl1",
                "Dnajc9","Nbeal1","Plod2","Ssh2","Trappc10","Ercc6","Fchsd2",
                "Gfi1","Ubn2","Vps54","Actr1b","Ccni","Cd2bp2","Tnfsf10",
                "Acot11","Atad2","Lgals9","Nup153","Gtpbp1")
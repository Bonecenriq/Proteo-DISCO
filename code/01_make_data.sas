/*==============================================================================
  STUDY       : SOMMA
  PROGRAM     : 01_make_data.sas

  PURPOSE     : Create the analytic datasets used for the proteomics entropy
                analyses.

  INPUTS      : SOMMA Parent and SOMMA Jr baseline datasets, including
                demographic, clinical, physical function, body composition,
                and proteomics-related variables.

  OUTPUTS     : DATALIB.PROT
                DATALIB.PROT_SUBSET
                DATALIB.PROT_SUBSET_MEN
                DATALIB.PROT_SUBSET_WOMEN
                DATALIB.PROTJR
                DATALIB.PROTJR_SUBSET
                DATALIB.PROTJR_SUBSET_MEN
                DATALIB.PROTJR_SUBSET_WOMEN

  DATA ACCESS : Individual-level SOMMA data are not included in this
                repository. Access to the underlying study data is required
                to reproduce the analyses.

  DEPENDENCIES:
                - SOMMA baseline datasets
                - SOMMA Jr baseline datasets
                - SOMMA SAS formats

  VERSION     : 1.0
  DATE        : 2025-11-24
  AUTHOR      : T. Mansfield

  NOTES       : Local file paths are specified separately from the analysis
                code and should be configured for the user's environment.
==============================================================================*/

/*==============================================================================
  LOCAL CONFIGURATION
  ------------------------------------------------------------------------------

  Update the paths below for the local environment in which this program is
  being run.

  These paths should point to the required SOMMA data, SOMMA Jr data, SAS
  formats, and external macro files.

  The paths below are intentionally not populated with institution-specific
  locations so that this program can be shared publicly.
==============================================================================*/

%LET PATH_BASELINE = /path/to/somma/baseline;
%LET PATH_BASEJR   = /path/to/somma/jr/baseline;
%LET PATH_BASEJR2  = /path/to/somma/jr/baseline;
%LET PATH_FORMATS  = /path/to/somma/formats;
%LET PATH_DATA     = /path/to/analysis/data;


/*------------------------------------------------------------------------------
  SAS OPTIONS
------------------------------------------------------------------------------*/

OPTIONS MPRINT NOFMTERR MERGENOBY=WARN;


/*------------------------------------------------------------------------------
  LIBRARIES
------------------------------------------------------------------------------*/

LIBNAME BASELINE "&PATH_BASELINE";
LIBNAME BASEJR   "&PATH_BASEJR";
LIBNAME BASEJR2  "&PATH_BASEJR2";
LIBNAME LIBRARY  "&PATH_FORMATS";
LIBNAME DATALIB  "&PATH_DATA";


/*==============================================================================
  SOMALOGIC ANALYTIC SUBSET
==============================================================================*/

* This file simply has the SOMMA IDs for those included in the SomaScan dataset *;
FILENAME CSV "&PATH_DATA/somalogic_subset_ids.csv"
    TERMSTR=CRLF;

PROC IMPORT DATAFILE=CSV
    OUT=WORK.SOMALOGIC_IDS
    DBMS=CSV
    REPLACE;
RUN;


*===================================================================================*;
* SOMMA PARENT *
*===================================================================================*;

DATA PROT; 
    MERGE BASELINE.VV000_2025_11_W1(KEEP=ID SITE
					/*** DEMOGRAPHICS ***/
						D1AGE2 EEFEMALE EERACE EERCWRMH 
						EDEDUC TURSMOKE TUDRPRWK
					/*** CARDIOVASCULAR RISK FACTORS ***/
						HWBMI HCWAIST HWHGT HWWGT
					/*** TABLE 1 VARIABLES ***/
						PACALWKALL MHSOMMAMI NF400MPACE NFSPPBE
						QLHEALTH GSGRPMAX MHHEALTH
					/*** SECONDARY OUTCOMES ***/
						DSCORR
					/*** CARDIOVASCULAR RISK FACTORS ***/
						MHMI MHCHF MHAFIB MHASTEN MHPVD MHSTRK MHBLKART 
						
					/*** OTHER ORGAN SPECIFIC DISEASES ***/
						MHRENAL MHKC MHCOPD MHLC MHOSTEO
						
						IN = INBASE)
		  BASELINE.AC000_2025_11_W1(KEEP = ID ACTOTCT24IM)
		  BASELINE.ACGG000_2025_12_W1(KEEP = ID ACSCFUAV)
		  BASELINE.FT000_2026_01_W1(KEEP = ID FT0V2FN)
		  BASELINE.MR000_2025_05_W1(KEEP = ID MRMAFATTOT MRMTFFM MRMLQFFM MRMRQFFM)
		  			/*** TABLE 1 VARIABLES ***/					
		  BASELINE.CR000_2025_07_W1(KEEP = ID CRMMASSSC CRHID3CR)
		  BASELINE.LE000_2024_11_W1(KEEP = ID LEPEAKPWR2 LEREPMAX LELEG)
		  BASELINE.TT000_2024_11_W1(KEEP = ID TTPKVO2 TTPKVO2U TTSSPKVO2)
		  BASELINE.RE000_2022_04_W1(KEEP = ID REMOXPHOS REP1TECH)
		  BASELINE.M1000_2025_05_W1(KEEP = ID M1ANMEDS)
		  SOMALOGIC_IDS(IN=INSOMALOGIC)
		  ; 
			
	BY ID;
	IF INBASE;
	
	/*** INSUBSET: OUR ANALYTIC SUBSET IS ALL OF THOSE WITH SOMALOGIC DATA ***/
	IF INSOMALOGIC = 1 THEN INSUBSET = 1;
	ELSE INSUBSET = 0;
	
	FORMAT INSUBSET YNDK.;
	LABEL INSUBSET = "IN ANALYTIC SUBSET?";
RUN;
	

/* 	Make analytic dataset */
DATA DATALIB.PROT;
    SET PROT;
    
    /* HEIGHT IN CM */
    HWHGTCM = HWHGT * 100;
    LABEL HWHGTCM = "Height in centimeters";
	FORMAT HWHGTCM SPECF.;
   	
   	
   	/* SELF-REPORTED HEALTH */
   IF MHHEALTH = 1 THEN DO; MHHEALTHCAT = 3; MHHEALTHCAT1 = 0; MHHEALTHCAT2 = 0; MHHEALTHCAT3 = 1; END;
   ELSE IF MHHEALTH = 2 THEN DO; MHHEALTHCAT = 2; MHHEALTHCAT1 = 0; MHHEALTHCAT2 = 1; MHHEALTHCAT3 = 0; END;
   ELSE IF MHHEALTH >= 3 THEN DO; MHHEALTHCAT = 1; MHHEALTHCAT1 = 1; MHHEALTHCAT2 = 0; MHHEALTHCAT3 = 0; END;
   ELSE DO; MHHEALTHCAT = MHHEALTH; MHHEALTHCAT1 = MHHEALTH; MHHEALTHCAT2 = MHHEALTH; MHHEALTHCAT3 = MHHEALTH; END;
   
   LABEL MHHEALTHCAT = "SELF-REPORTED HEALTH CATEGORY: GOOD/FAIR, VERY GOOD, EXCELLENT"
   		MHHEALTHCAT1 = "SELF-REPORTED HEALTH GOOD/FAIR"
   		MHHEALTHCAT2 = "SELF-REPORTED HEALTH VERY GOOD"
   		MHHEALTHCAT3 = "SELF-REPORTED HEALTH EXCELLENT";
   		
   	/* 	MULTIMORBIDITY */
	IF MHSOMMAMI <= .Z THEN MULTMORBCAT=MHSOMMAMI;
	ELSE MULTMORBCAT = MIN(MHSOMMAMI, 2);
	FORMAT MULTMORBCAT SPECF.;
	LABEL MULTMORBCAT="Multimorbidity category (0, 1, 2+)";
	
	/*** CATEGORICAL EDUCATION VARIABLE ***/
	IF EDEDUC>.Z THEN DO;
		EDUC1=0; EDUC2=0; EDUC3=0;
		
		*HIGH SCHOOL OR LESS OR OTHER;
		IF 0<=EDEDUC<=2 OR EDEDUC=7 THEN DO; EDUCCAT=1; EDUC1=1; END;
		IF EDEDUC=3 THEN DO; EDUCCAT=2; EDUC2=1; END; *SOME COLLEGE;
		IF EDEDUC=4 OR EDEDUC=5 THEN DO; EDUCCAT=3; EDUC3=1; END; *COLLEGE GRADUATE OR MORE;
	END;
	ELSE DO;
		EDUC1=EDEDUC; EDUC2=EDEDUC; EDUC3=EDEDUC;
		EDUCCAT=EDEDUC;
	END;
	
	LABEL EDUCCAT="Education level 3 cats"
		EDUC1="High school or less or other "
		EDUC2="Some college"
		EDUC3="College graduate or more";
	FORMAT EDUC1 EDUC2 EDUC3 YNDK. EDUCCAT MYEDUCF.;

	
	
	/*** ALCOHOL AND SMOKING ***/
	IF TURSMOKE = 0 THEN SMOKECAT = 0;
	IF TURSMOKE = 1 OR TURSMOKE = 2 THEN SMOKECAT = 1;
	LABEL SMOKECAT="Current/Past Smoker";
	FORMAT SMOKECAT YNDK.;
	
	IF TUDRPRWK>.Z THEN DO;
		ALC1=0; ALC2=0; ALC3=0;
		IF TUDRPRWK = 0 THEN DO; 	 ALCCAT=1; ALC1=1; END; *NONDRINKER;
		IF 0 < TUDRPRWK < 7 THEN DO; ALCCAT=2; ALC2=1; END; *LIGHT DRINKER;
		IF 7 <= TUDRPRWK THEN DO; 	 ALCCAT=3; ALC3=1; END; *HEAVY DRINKER;
	END;
	ELSE DO;
		ALC1=TUDRPRWK; ALC2=TUDRPRWK; ALC3=TUDRPRWK; ALCCAT=TUDRPRWK;
	END;
	
	LABEL ALCCAT="Alcohol Consumption Category (1-3)" 
		  ALC1 = "Nondrinker"
		  ALC2 = "Light Drinker (1-<7 drinks per week)"
		  ALC3 = "Heavy drinker (>=7 drinks per week)";
	FORMAT ALC1 ALC2 ALC3 YNDK.;
	
	/*** Steps per day in 1000s ***/
	ACSCFUAVK = ACSCFUAV / 1000;
	LABEL ACSCFUAVK = "Mean step count (thousand steps/day)";
	FORMAT ACSCFUAVK SPECF.;
	
	/* EXCLUDE PARTICIPANTS WITH UNUSUAL D3CR VALUES FROM MUSCLE MASS */
	IF CRHID3CR = 1 THEN CRMMASSSC = .W;
	
	IF MHMI = 1 OR MHCHF = 1 OR MHSTRK = 1 THEN MHCVD = 1;
	ELSE IF MHMI = 0 AND MHCHF = 0 AND MHSTRK = 0 THEN MHCVD = 0;
	ELSE MHCVD = .A;
	IF MHKC = .M THEN MHKC = 0;
	IF MHLC = .M THEN MHLC = 0;
RUN;

PROC FREQ DATA=DATALIB.PROT;
	TABLES MHHEALTH*MHHEALTHCAT*MHHEALTHCAT1*MHHEALTHCAT2*MHHEALTHCAT3 
	MHMI*MHCHF*MHSTRK*MHCVD
	MHMI MHCHF MHAFIB MHASTEN MHPVD MHSTRK MHBLKART 
	MHRENAL MHKC MHCOPD MHLC MHOSTEO/ LIST MISSPRINT;


/* 	Save subsets of the data*/
DATA DATALIB.PROT_SUBSET;
	SET DATALIB.PROT;
	IF INSUBSET = 1;
RUN;

DATA DATALIB.PROT_SUBSET_MEN;
	SET DATALIB.PROT_SUBSET;
	IF EEFEMALE = 0;
RUN;

DATA DATALIB.PROT_SUBSET_WOMEN;
	SET DATALIB.PROT_SUBSET;
	IF EEFEMALE = 1;
RUN;


*===================================================================================*;
* SOMMA JR *
*===================================================================================*;

DATA PROTJR; 
    MERGE BASEJR.JRVV000_2025_11(KEEP=ID SITE
					/*** DEMOGRAPHICS ***/
						D1AGE2 EEFEMALE EERACE EERCWRMH 
						EDEDUC TURSMOKE TUDRPRWK
					/*** CARDIOVASCULAR RISK FACTORS ***/
						HWBMI HWHGT HWWGT
					/*** TABLE 1 VARIABLES ***/
						MHSOMMAMI NF400MPACE NFSPPBE
						QLHEALTH GSGRPMAX MHHEALTH DSCORR

						IN = INBASE)
		  BASEJR2.JRACGG000_2026_05(KEEP = ID ACSCFUAVM RENAME=(ACSCFUAVM=ACSCFUAV))
		  BASEJR.JRMR000_2024_05(KEEP = ID MRMAFATTOT MRMTFFM MRMLQFFM MRMRQFFM)
		  			/*** TABLE 1 VARIABLES ***/					
		  BASEJR.JRCR000_2025_07(KEEP = ID CRMMASSSC CRHID3CR)
		  BASEJR.JRLE000_2024_11(KEEP = ID LEPEAKPWR2 LEREPMAX LELEG)
		  BASEJR.JRTT000_2024_11(KEEP = ID TTPKVO2 TTPKVO2U)
		  BASEJR.JRRE000_2022_11(KEEP = ID REMOXPHOS)
		  SOMALOGIC_IDS(IN=INSOMALOGIC)
		  ; 
			
	BY ID;
	IF INBASE;
	
	/*** OUTCOME FLAGS ***/
	IF INSOMALOGIC = 1 THEN INSUBSET = 1;
	ELSE INSUBSET = 0;
RUN;
	

/* 	Make analytic dataset */
DATA DATALIB.PROTJR;
    SET PROTJR;
    
    /* HEIGHT IN CM */
    HWHGTCM = HWHGT * 100;
   	
   	
   	/* SELF-REPORTED HEALTH */
   IF MHHEALTH = 1 THEN DO; MHHEALTHCAT = 3; MHHEALTHCAT1 = 0; MHHEALTHCAT2 = 0; MHHEALTHCAT3 = 1; END;
   ELSE IF MHHEALTH = 2 THEN DO; MHHEALTHCAT = 2; MHHEALTHCAT1 = 0; MHHEALTHCAT2 = 1; MHHEALTHCAT3 = 0; END;
   ELSE IF MHHEALTH >= 3 THEN DO; MHHEALTHCAT = 1; MHHEALTHCAT1 = 1; MHHEALTHCAT2 = 0; MHHEALTHCAT3 = 0; END;
   ELSE DO; MHHEALTHCAT = MHHEALTH; MHHEALTHCAT1 = MHHEALTH; MHHEALTHCAT2 = MHHEALTH; MHHEALTHCAT3 = MHHEALTH; END;
   
   LABEL MHHEALTHCAT = "SELF-REPORTED HEALTH CATEGORY: GOOD/FAIR, VERY GOOD, EXCELLENT"
   		MHHEALTHCAT1 = "SELF-REPORTED HEALTH GOOD/FAIR"
   		MHHEALTHCAT2 = "SELF-REPORTED HEALTH VERY GOOD"
   		MHHEALTHCAT3 = "SELF-REPORTED HEALTH EXCELLENT";
   		
   	/* 	MULTIMORBIDITY */
	IF MHSOMMAMI <= .Z THEN MULTMORBCAT=MHSOMMAMI;
	ELSE MULTMORBCAT = MIN(MHSOMMAMI, 2);
	FORMAT MULTMORBCAT SPECF.;
	LABEL MULTMORBCAT="Multimorbidity category (0, 1, 2+)";
	
	/*** CATEGORICAL EDUCATION VARIABLE ***/
	IF EDEDUC>.Z THEN DO;
		EDUC1=0; EDUC2=0; EDUC3=0;
		
		*HIGH SCHOOL OR LESS OR OTHER;
		IF 0<=EDEDUC<=2 OR EDEDUC=7 THEN DO; EDUCCAT=1; EDUC1=1; END;
		IF EDEDUC=3 THEN DO; EDUCCAT=2; EDUC2=1; END; *SOME COLLEGE;
		IF EDEDUC=4 OR EDEDUC=5 THEN DO; EDUCCAT=3; EDUC3=1; END; *COLLEGE GRADUATE OR MORE;
	END;
	ELSE DO;
		EDUC1=EDEDUC; EDUC2=EDEDUC; EDUC3=EDEDUC;
		EDUCCAT=EDEDUC;
	END;
	
	LABEL EDUCCAT="Education level 3 cats"
		EDUC1="High school or less or other "
		EDUC2="Some college"
		EDUC3="College graduate or more";
	FORMAT EDUC1 EDUC2 EDUC3 YNDK. EDUCCAT MYEDUCF.;

	
	
	/*** ALCOHOL AND SMOKING ***/
	IF TURSMOKE = 0 THEN SMOKECAT = 0;
	IF TURSMOKE = 1 OR TURSMOKE = 2 THEN SMOKECAT = 1;
	LABEL SMOKECAT="Current/Past Smoker";
	FORMAT SMOKECAT YNDK.;
	
	IF TUDRPRWK>.Z THEN DO;
		ALC1=0; ALC2=0; ALC3=0;
		IF TUDRPRWK = 0 THEN DO; 	 ALCCAT=1; ALC1=1; END; *NONDRINKER;
		IF 0 < TUDRPRWK < 7 THEN DO; ALCCAT=2; ALC2=1; END; *LIGHT DRINKER;
		IF 7 <= TUDRPRWK THEN DO; 	 ALCCAT=3; ALC3=1; END; *HEAVY DRINKER;
	END;
	ELSE DO;
		ALC1=TUDRPRWK; ALC2=TUDRPRWK; ALC3=TUDRPRWK; ALCCAT=TUDRPRWK;
	END;
	
	LABEL ALCCAT="Alcohol Consumption Category (1-3)" 
		  ALC1 = "Nondrinker"
		  ALC2 = "Light Drinker (1-<7 drinks per week)"
		  ALC3 = "Heavy drinker (>=7 drinks per week";
	FORMAT ALC1 ALC2 ALC3 YNDK.;
	
	/* EXCLUDE PARTICIPANTS WITH UNUSUAL D3CR VALUES FROM MUSCLE MASS */
	IF CRHID3CR = 1 THEN CRMMASSSC = .W;
RUN;

PROC FREQ DATA=DATALIB.PROTJR;
	TABLES MHHEALTH*MHHEALTHCAT*MHHEALTHCAT1*MHHEALTHCAT2*MHHEALTHCAT3 / LIST MISSPRINT;

/* 	Save subsets of the data*/

DATA DATALIB.PROTJR_SUBSET;
	SET DATALIB.PROTJR;
	IF INSUBSET = 1;
RUN;

DATA DATALIB.PROTJR_SUBSET_MEN;
	SET DATALIB.PROTJR_SUBSET;
	IF EEFEMALE = 0;
RUN;

DATA DATALIB.PROTJR_SUBSET_WOMEN;
	SET DATALIB.PROTJR_SUBSET;
	IF EEFEMALE = 1;
RUN;
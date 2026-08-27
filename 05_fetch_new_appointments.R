# -------------------------------------------------------------------------
# 1. Setup and load configurations
# -------------------------------------------------------------------------
source("00_libraries_and_utils.R")
# Load configuration parameters
conf <- config::get()

new_appointments <- local({
  con <- dbConnect(
    odbc(),
    Driver = "ODBC Driver 17 for SQL Server",
    Server = "BIAtlasdev\\BIAtlasdev",
    Database = "Data_Analytics_Power_BI",
    Trusted_Connection = "Yes",
    TrustServerCertificate = "yes"
  )
  
  dataset <- dbGetQuery(
    con,
    "


-- 18/08/26 KL Urgency changed to referral - Elise Sapsford

--WITH cte_distance AS
--	(
--SELECT [Postcode_8_chars] AS sm1
--      ,[Latitude_1m] AS sm_lat
--      ,[Longitude_1m] AS sm_long
--  FROM [UK_Health_Dimensions].[ODS].[Postcode_Grid_Refs_Eng_Wal_Sco_And_NI_SCD]
--  WHERE Postcode_8_chars IN ('BS10 5NB')
--  AND Effective_To IS NULL
--	)

DECLARE @RunDate DATE = GETDATE();
SET DATEFIRST 1;


WITH ClinicSites AS (
    SELECT *
    FROM (
        VALUES
            ('Southmead Hospital', 'BS10 5NB'),
            ('COSSHAM HOSPITAL', 'BS15 1LF'),
            ('YATE HEALTH CENTRE', 'BS37 4BQ'),
            ('DOWNEND CLINIC', 'BS16 5SG'),
            ('HORFIELD HEALTH CENTRE', 'BS7 0XW'),
            ('WESTON GENERAL HOSPITAL', 'BS23 4TQ'),
            ('CADBURY HEATH HEALTH CENTRE', 'BS30 8EL'),
            ('Non NHS Location', 'BS10 5NB'),
            ('THORNSBURY CLDT', 'BS35 2AB'),
            ('SOUTHMEAD HEALTH CENTRE', 'BS10 5BP'),
            ('BATH OUTPATIENT CLINIC', 'BA1 3NG'),
            ('FISHPONDS HEALTH CENTRE', 'BS16 3TD'),
            ('FRENCHAY HOSPITAL', 'BS16 1LE'),
            ('GLOUCESTER MEDICAL CENTRE', 'BS7 8PS'),
            ('WELLSPRING HEALTHY LIVING CENTRE', 'BS5 9QY'),
            ('BRISTOL ROYAL INFIRMARY', 'BS2 8HW'),
            ('CLEVEDON HEALTH CENTRE', 'BS21 6DG'),
            ('SHIREHAMPTON HEALTH CENTRE', 'BS11 9SB'),
            ('JOHN MILTON CLINIC', 'BS10 7DP'),
            ('ST GEORGE HEALTH CENTRE', 'BS5 7PH'),
            ('SOUTH BRISTOL RENAL UNIT', 'BS4 1WH'),
            ('CLEVEDON HOSPITAL', 'BS21 7HN'),
            ('GLOUCESTER OUTPATIENT CLINIC', 'GL1 3NN'),
            ('CHIPPENHAM COMMUNITY HEALTH', 'SN15 2AJ'),
            ('BATH RENAL OUT POST', 'BA1 3NG'),
            ('WESTON RENAL OUT POST', 'BS23 4TQ'),
            ('COAST RESOURCE CENTRE', 'BS20 7AW'),
            ('FILTON CLINIC', 'BS34 7BQ')
    ) AS v (ClinicSite, ClinicPostcode)
)
, cte_distance AS (

    SELECT
   [Postcode_8_chars] AS sm1
      ,[Latitude_1m] AS sm_lat
      ,[Longitude_1m] AS sm_long
	  ,ClinicSite
	  ,ClinicPostcode

    FROM ClinicSites sp

    LEFT JOIN [UK_Health_Dimensions].[ODS].[Postcode_Grid_Refs_Eng_Wal_Sco_And_NI_SCD] pc
        ON Postcode_single_space_e_Gif  = sp.ClinicPostcode   AND Effective_To IS NULL

)



SELECT DISTINCT
    CAST(a.AppointmentDTTM AS DATE)                  AS appt_dttm,
    a.PATIENT_ID                       AS dim_patient_id,
a.LocalTreatmentFunctionCode       AS local_spec_code,
    a.TreatmentFunctionCode            AS national_spec_code,
	a.AppointmentTypeCode              AS appointment_type_code,
    a.AppointmentType                  AS appointment_type,
a.ConsultationMediaCode            AS consultation_media_code,
    a.ConsultationMedia                AS consultation_media,
a.AttendedStatusCode               AS attended_status_code,
    a.AttendedStatus                   AS attended_status,


  6371 * ACOS(COS(RADIANS(r.sm_lat)) * COS(RADIANS(ps.[Latitude_1m])) *
				COS(RADIANS(ps.[Longitude_1m] - r.sm_long)) +
				SIN(RADIANS(r.sm_lat)) * SIN(RADIANS(ps.[Latitude_1m]))
                ) AS distance_km,

    CASE WHEN P.PatientPostcode = 'ZZ99 3VZ' THEN 1 ELSE 0 END AS nfa_ind   ,
	CASE
	WHEN a.AgeAtAppointment  <18 THEN '<18'
	WHEN a.AgeAtAppointment  < 25 THEN '18-25'
	WHEN AgeAtAppointment < 35 THEN '25-34'
	WHEN AgeAtAppointment < 45 THEN '35-44'
	WHEN AgeAtAppointment < 55 THEN '45-54'
	WHEN AgeAtAppointment < 65 THEN '55-64'
	WHEN AgeAtAppointment < 75 THEN '65-74'
	WHEN AgeAtAppointment < 85 THEN '75-84'
	WHEN AgeAtAppointment < 95 THEN '85-94'
	ELSE '95+' END AS age_group,
    a.AgeAtAppointment                 AS age_at_appointment,
 p.EthnicityCode                    AS ethnicity_code,
    p.Ethnicity                        AS ethnicity,
    p.Index_Multiple_Deprivation_Decile AS index_multiple_deprivation_decile,


ISNULL(Alerts.a_ld, 0)                        AS a_ld,
ISNULL(Alerts.a_autism, 0)                   AS a_autism,
ISNULL(Alerts.a_interpreter_req_bsl, 0)      AS a_interpreter_req_bsl,
ISNULL(Alerts.a_interpreter_req_lang, 0)     AS a_interpreter_req_lang,
0                                            AS a_balance,
ISNULL(Alerts.a_cognitive_impairment, 0)     AS a_cognitive_impairment,
ISNULL(Alerts.a_mobility_restriction, 0)     AS a_mobility_restriction,
ISNULL(Alerts.a_hear_vis_impaired, 0)        AS a_hear_vis_impaired,
ISNULL(Alerts.a_dementia, 0)                 AS a_dementia,
0 AS a_depression,
0 AS a_downs_syndrome,
0 AS a_long_standing_condition,
0 AS a_makaton,
0 AS a_mild_cognitive_impairment,
0 AS a_memory_impairment,
0 AS a_mood_disorder,
0 AS a_other_disability,
0 AS a_psychosis,
0 AS a_severe_anxiety,
0 AS a_wheelchair_user,


    p.GenderCode                       AS gender_code,
    p.Gender                           AS gender,

    a.RegisteredGPPracticeCode         AS registered_gp_practice_code,
    a.RegisteredGPPractice             AS registered_gp_practice,

    a.SiteCode                         AS site_code,

    ISNULL(HIST.dna_count, 0) AS prev_dna_LY,

    DATEPART(HOUR, a.AppointmentDTTM)  AS appt_hour,
    DATENAME(WEEKDAY, a.AppointmentDTTM) AS appt_dow,
    DATEADD(month, DATEDIFF(month, 0, a.AppointmentDTTM), 0)          AS appt_month,
    CASE
        WHEN DATENAME(WEEKDAY, a.AppointmentDTTM) IN ('Saturday', 'Sunday')
        THEN 1 ELSE 0
    END                                AS appt_wknd_ind,

    a.Urgency                          AS referral_urgency,
    DATEDIFF(DAY,a.AppointmentMadeDTTM,a.AppointmentDTTM)    AS lead_time_days,       -- derive if needed

    a.ClinicCode                       AS clinic_code,
    a.ClinicLocation                   AS clinic_location
	--,  CASE
 --          WHEN RAND(CHECKSUM(NEWID())) < 0.5 THEN 'Training'
 --          ELSE 'Testing' END AS test_train

,a.EXTERNAL_ID
,HomePhoneNumber
,MobilePhoneNumber
,ClinicName	,ClinicSessionCode	,ClinicSession
,p.Surname
,a.AppointmentDTTM AS [Appt Date/Time]
,LocalTreatmentFunction
,Site
,a.ClinicSite
,a.NHSNumber
,a.PASID
,p.DateOfBirth
,  NULL AS lsoa_code -- IMD.[LSOA code (2011)] AS lsoa_code
        ,CASE
            WHEN CAST(a.AppointmentDTTM AS DATE) > CAST(GETDATE() AS DATE) THEN  'Testing'
            ELSE 'Training'
        END AS test_train  -- testing = future, training = past
		,[ReferralUrgencyCode] as UrgencyCode
		,[ReferralUrgency] as Urgency
		,OP_APPT_ID
		,trim(left(P.PatientPostcode,4)) as PostalDistrict
		,case when  clinicSession in
(
'BBS F2F-Gastro Surgery Upper GI -AdHoc',
'BBS F2F-Gastro Surgery UpperGI Surg REG -MTWHF E1W',
'BBS F2F-Gastro Surgery Upper Gastro CNS -M E1W',
'BBS F2F - Gastro Surgery Upper Gastro CNS - M E1W', -- looks same except for spaces
'BBS F2F-Gastro Surgery UpperGI Surg REG -MTWHF E1W',
'BBS F2F-Gastro Surgery Upper GI CNS -T E1W'
) then 'Yes' else 'No' end as IncludeClinics

,CASE TreatmentFunctionCode
	WHEN '101' THEN 'Urology'
	WHEN '301' THEN 'Upper GI'
	WHEN '502' THEN 'Gynaecology'
	WHEN '103' THEN 'Breast'
	ELSE 'Unknown'
END AS TumourSite




FROM nhs_trust_careflow_reporting.CFL.tbl_OP_APPOINTMENTS A
left join [nhs_trust_careflow_reporting].[CFL].[tbl_OP_REFERRALS] ref on a.[REFERRAL_ID] = ref.[REFERRAL_ID]

LEFT JOIN nhs_trust_careflow_reporting.CFL.tbl_PATIENTS p ON A.PATIENT_ID = p.PATIENT_ID

LEFT JOIN [UK_Health_Dimensions].[ODS].[Postcode_Grid_Refs_Eng_Wal_Sco_And_NI_SCD] AS LSOA WITH (NOLOCK) ON REPLACE(UPPER(LSOA.Postcode_8_chars),' ','') = REPLACE(UPPER(P.PatientPostcode),' ','')
--LEFT JOIN [nhs_trust_careflow_semantic_layer].[Reference].[tbl_IMD2019_Index_of_Multiple_Deprivation] AS IMD WITH (NOLOCK) ON IMD.[LSOA code (2011)] = LSOA.LSOA

LEFT JOIN [UK_Health_Dimensions].[ODS].[Postcode_Grid_Refs_Eng_Wal_Sco_And_NI_SCD] ps WITH (NOLOCK) ON REPLACE(UPPER(ps.Postcode_8_chars),' ','') = REPLACE(UPPER(P.PatientPostcode),' ','') AND  ps.Effective_To IS NULL

left join cte_distance r on a.clinicsite = r.ClinicSite

--CROSS JOIN cte_distance r

LEFT JOIN (
    SELECT
        PATIENT_ID,
        MAX(CASE WHEN Alert LIKE 'Learning Disability' THEN 1 ELSE 0 END) AS a_ld,
        MAX(CASE WHEN Alert LIKE 'Autistic Spectrum Disorder' THEN 1 ELSE 0 END) AS a_autism,
        MAX(CASE WHEN Alert LIKE 'Communication Support Used - Sign Language' THEN 1 ELSE 0 END) AS a_interpreter_req_bsl,
        MAX(CASE WHEN Alert LIKE 'Interpreter required' THEN 1 ELSE 0 END) AS a_interpreter_req_lang,
        MAX(CASE WHEN Alert = 'Cognitive Impairment' THEN 1 ELSE 0 END) AS a_cognitive_impairment,
        MAX(CASE WHEN Alert = 'Mobility Impairment' THEN 1 ELSE 0 END) AS a_mobility_restriction,
        MAX(CASE WHEN Alert = 'Hearing Impaired' THEN 1 ELSE 0 END) AS a_hear_vis_impaired,
        MAX(CASE WHEN Alert = 'Cognitive Impairment - Dementia' THEN 1 ELSE 0 END) AS a_dementia
    FROM [nhs_trust_careflow_reporting].[CFL].[tbl_PATIENT_ALERTS]
    WHERE AlertEndDate IS NULL
    GROUP BY PATIENT_ID
) Alerts
    ON p.PATIENT_ID = Alerts.PATIENT_ID

			OUTER APPLY (
					    SELECT COUNT(*) AS dna_count
	    FROM  nhs_trust_careflow_reporting.CFL.tbl_OP_APPOINTMENTS AS HIST_OP WITH (NOLOCK)
	    WHERE HIST_OP.PATIENT_ID = a.PATIENT_ID
	      AND HIST_OP.AppointmentDTTM < a.AppointmentDTTM  -- Only count past appointments
	      -- Use IN here to handle multiple IDs mapped to NHS_ID '3'
	      AND HIST_OP.AttendedStatus IN ('DNA (Not Specified)')
	      AND HIST_OP.AppointmentDTTM >= DATEADD(YEAR, -1, a.AppointmentDTTM)
		  ) AS hist

			WHERE 1=1   
			-- w/c 2 weeks time
			and CAST(a.AppointmentDTTM AS DATE) >= DATEADD(WEEK, 2, DATEADD(DAY, 1 - DATEPART(WEEKDAY, @RunDate), CAST(@RunDate AS DATE)))  -- w/c 2 weeks time
			AND CAST(a.AppointmentDTTM AS DATE) < DATEADD(DAY, 7,DATEADD(WEEK, 2, DATEADD(DAY, 1 - DATEPART(WEEKDAY, @RunDate), CAST(@RunDate AS DATE))))			
			
			--CAST(a.AppointmentDTTM AS DATE)>= cast(getdate() as date)
			--and  CAST(a.AppointmentDTTM AS DATE)<= cast(getdate()+21 as date)
			AND a.TreatmentFunctionCode IN ('101','301','502','103')
			and [ReferralUrgencyCode] = '3'
			and AttendedStatusCode = 0"
  )
dbDisconnect(con)
dataset
})

write.csv(new_appointments, conf$new_appointment_path)

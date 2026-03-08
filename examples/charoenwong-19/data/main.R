#####################################################################
# Date: April 23, 2019
# Description:
# Code to run the main analyses in 
# Charoenwong, Kwan, Umar. "Does Regulatory Jurisdiction Affect the 
#      Quality of Investment-Adviser Regulation?". American Economic 
#      Review. Forthcoming.
#
# Note that we do not include code for the additional analyses found 
# in the appendix. However, they are all available upon request.
#
# Please get in touch at any of the following emails
# if you have any questions:
# 1. Alan Kwan <apkwan@hku.hk> | <alanpaulkwan@gmail.com>
# 2. Ben Charoenwong <bizbgc@nus.edu.sg> | <bcharoen@chicagobooth.edu>
# 3. Tarik Umar <tarik.umar@rice.edu> | <tarik.umar@chicagobooth.edu>
#
######################################################################
# Date: August 18, 2024
# Description: Modified Step #1 and first few lines in Step #2
# to control for transitions appropriately
# Modified by Young Ahn <youngahn@sas.upenn.edu>
######################################################################

############## Step #1 : read helper funcitons functions (changed by Young Ahn)  ###############
  # Note: adjust the working directory reference folder as necessary
folder <- paste0(getwd(),"/")
out.this='test_run/'  # Where outputs will go
dir.create(out.this,showWarnings = FALSE)
source('code_for_distribution/functions_read_in.R') # Reads in all functions

############## Step #2 : read in data         ###############

for(item in list.files(paste(folder,'data_for_distribution/',sep=''),full.names = TRUE,pattern='fst')){
   print(paste(item,'---',Sys.time()))
   # Read in, automatically assigned name
      #strips out the prefix 
    auto_name = strsplit(item,split='/') %>% sapply(dplyr::last) %>% gsub('\\.fst','',.)
    
    read_fst(item,as.data.table = TRUE) %>% assign(x=auto_name,value=,.,envir=.GlobalEnv) 
}

bdCRDs=readRDS('data_for_distribution/bd_crds.rds') # those whose ADV indicaters a broker dealer service
states=read_rds('data_for_distribution/state_crds.rds')
individual_panel[,treated:=treated1][,post:=year>=2012]
firm_panel[,treated:=treated1][,post:=year>=2012]

############## Step #3 : flip on switches to run diff analyses   ###############


# Note: All switches with the "_extra" suffix are not shown in the main paper tables. 
# but are provided as additional robustness checks for our results.
# Set to be TRUE to show additional robustness checks and specifications.

runMain=TRUE

runAccusationTypes=TRUE
runAccusationTypes_extra=FALSE
  
runBudgetSalary=TRUE
runBudgetSalary_extra=FALSE

runDistance=TRUE
runDistance_extra=FALSE

runRegAction=TRUE
runRegAction_extra=FALSE

runIndivMarketPower=TRUE
runIndivMarketPower_extra=FALSE

runRepeatComplaints=TRUE
runRepeatComplaints_extra=TRUE

runDualRegistered=TRUE
runDualRegistered_extra=FALSE

runPTgraphs=TRUE
runRobustness=TRUE

############### Step #4 : run switches    ##################
#   Note: the code was refactored for documentation purposes only

# Table 2: Main Result----
if(runMain){
  
  # firm - level
  {
    l=list()
    
    firmsubset=firm_panel[year %in% seq(2009,2014),]
    
    felm(data=firmsubset,I(100*winsor(n_complaints/N,h=0.98)) ~ post*treated+log(N)+I(log(N)^2)+
           I(log(N)^3)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    
    felm(data=firmsubset,I(log1(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)
         |firmcrd+year|0|state.final) %>% summary2
    felm(data=firmsubset,I(log1(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)
         |firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    
    felm(data=firmsubset,I(log1(100*n_complaints/N)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)
         |firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    
    felm(data=firmsubset,I(sign2(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    # dont record
    felm(data=firmsubset,I(sign2(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|
           firmcrd+I(paste(state.final,year))|0|state.final) %>% summary
    felm(data=firmsubset[year %in% seq(2008,2015),],I(log1(n_complaints)) ~ post*treated+log(N)+
           I(log(N)^2)+I(log(N)^3)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    felm(data=firmsubset[year %in% seq(2010,2013),],I(log1(n_complaints)) ~ post*treated+log(N)+
           I(log(N)^2)+I(log(N)^3)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary
    
    stargazer(l,type='text',omit.stat='ser',keep='treated',
              add.lines=list(c('state fe',grepfe(l,'state')),
                             c('firm fe',grepfe(l,'firm')),
                             c('period','9/14','9/14','9/14','9/14','9/14','8/15')
              ));
    print("Column (1) not included in Table 2 for space. Columns (2) through (6) match Table 2 Panel A Columns (1) through (5).");
    
    rm(firmsubset) ; gc();
    
  } 
  # individual level
  {
    
    indivsubset=copy(individual_panel)
    setkey(indivsubset,crd,year)
    indivsubset[,complaint.lag.ever:=cumsum(coalesce(complaint.lag,0)) %>% sign,by=list(crd)]
    
    # lag complaint, state year, branch-post, branch-post * state-year, size-class 
    l=list()
    felm(data=indivsubset[year %in% seq(2009,2014),],I(sign2(complaint.sign)) ~ treated*post|firmcrd+year|0|state.final) %>% summary2
    felm(data=indivsubset[year %in% seq(2009,2014),],I(sign2(complaint.sign)) ~ treated*post|firmcrd+I(paste(state.final,post))|0|state.final) %>% summary2
    felm(data=indivsubset[year %in% seq(2009,2014),],I(sign2(complaint.sign)) ~ treated*post|firmcrd+I(paste(branchclean,post))+year|0|state.final) %>% summary2
    felm(data=indivsubset[year %in% seq(2009,2014),],I(sign2(complaint.sign)) ~ treated*post|firmcrd+I(paste(branchclean,post))+I(paste(state.final,year))|0|state.final) %>% summary2
    felm(data=indivsubset[year %in% seq(2009,2014),],I(sign2(complaint.sign)) ~ treated*post
         |firmcrd+crd+year|0|state.final) %>% summary2
    
        a=list(c('state fe',grepfe(l,'state')),
           c('firm fe',grepfe(l,'firm')),
           c('branch fe',grepfe(l,'branch')),
           c('indiv fe',grepfe(l,'indivcrd'))
    )
    stargazer(l,type='text',omit.stat='ser',
              add.lines=a,keep='treated')   
  }
}

# Table 3: Robustness----
if(runRobustness){
  
  firmsubset=firm_panel[year %in% seq(2008,2015),] 
  
  l=list()
  
  # at least 3 people
  felm(data=firmsubset[year %in% setdiff(seq(2009,2014),2012) & N>=1,],I(log1(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|
         firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  
  felm(data=firmsubset[year %in% seq(2009,2014) & N>=3,],I(log1(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|
         firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  
  felm(data=firmsubset[year %in% seq(2009,2014) & N<=50,],I(log1(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|
         firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  
  # no new york /wy
  felm(data=firmsubset[year %in% seq(2009,2014) & !(state.final %in% c('NY','WY')),],
       I(log1(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  # also no CA
  felm(data=firmsubset[year %in% seq(2009,2014) & !(state.final %in% c('NY','WY','CA')),],I(log1(n_complaints)) ~ 
         post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  
  
  # smaller
  felm(data=firmsubset[year %in% seq(2009,2014) & aum2<=(3*10^8),],I(log1(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|
         firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  felm(data=firmsubset[year %in% seq(2009,2014) & aum2<=(1*10^8),],I(log1(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|
         firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  
  felm(data=firmsubset[year %in% seq(2009,2014) & !(state.final %in% c('NY','WY','CA')),],
       I(sign2(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  felm(data=firmsubset[year %in% seq(2008,2015) & N<=50 & year!=2012,],I(sign2(n_complaints)) ~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|
         firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  
  stargazer(l,type='text',omit.stat='ser',no.space=TRUE,keep='treated',
            column.labels = c('no 2012','atl3','lt50','no ca/wy/ny','no ca/wy/ny','sub300','sub100','no2012,2008-2015,lt50'),
            notes = c("Columns (8) and (9) in this Table is not reported in Table 3.",
                      "Table 3 Column (8) in the paper repeats the Main results for comparison purposes."))
  
}

# Table 4: Complaint Type----
if(runAccusationTypes){
  x6=accusations_individual_year %>%
    merge(individual_panel,by=c('crd','year')) %>% 
    mutate(accuse_fraud = accuse_fraid) %>% setDT;
    # fixing previous typo.
  
  print("IAR (Individual)-level Results.");
  
  l=list()  
  
  
  x6[,fe2:=paste(branchclean,post)]
  felm(data=x6[year %in% seq(2009,2014),], sign2(accuse_fiduciary) ~ post*treated|firmcrd+fe2+year|0|state.final) %>% summary2
  felm(data=x6[year %in% seq(2009,2014),], sign2(accuse_fraud) ~ post*treated|firmcrd+fe2+year|0|state.final) %>% summary2
  felm(data=x6[year %in% seq(2009,2014),], sign2(accuse_churning) ~ post*treated|firmcrd+fe2+year|0|state.final) %>% summary2
  felm(data=x6[year %in% seq(2009,2014),], sign2(accuse_misrep) ~ post*treated|firmcrd+fe2+year|0|state.final) %>% summary2
  felm(data=x6[year %in% seq(2009,2014),], sign2(accuse_suit) ~ post*treated|firmcrd+fe2+year|0|state.final) %>% summary2
  
  
  # felm(data=x6[year %in% seq(2009,2014),], sign2(accuse_n5000) ~ post*treated|firmcrd+fe2+year|0|state.final) %>% summary2
    # Note: uncomment to show results. Significant but not included in Table 4 in the paper for space.
  felm(data=x6[year %in% seq(2009,2014),], sign2(accuse_n10000) ~ post*treated|firmcrd+fe2+year|0|state.final) %>% summary2
  # felm(data=x6[year %in% seq(2009,2014),], sign2(accuse_n50000) ~ post*treated|firmcrd+fe2+year|0|state.final) %>% summary2
    # Note: uncomment to show results. Significant but not included in Table 4 in the paper for space.
  felm(data=x6[year %in% seq(2009,2014),], sign2(accuse_n100000) ~ post*treated|firmcrd+fe2+year|0|state.final) %>% summary2
  
  l %>% stargazer(type='text',omit.stat='ser',
                  add.lines=list(c('firm fe?',grepfe(l,'firm')),
                                 c('year fe?',grepfe(l,'year')),
                                 c('branchpost fe?',grepfe(l,'branch')
                                 )),keep=':') 
  
  # firm level version
  if (runAccusationTypes_extra) {
    print("Firm-level Results.");
    print("Note: This is not shown in the paper, but shown here for robustness.");
    firm_panel2=x6[,list(
      N=sum(!is.na(complaint.sign)),
      state.final=state.final[1],
      n_complaints=sum(complaint.sign,na.rm=TRUE),
      n_complaints_total=sum(complaint.sign,na.rm=TRUE),
      n_reg=sum(complaint.reg,na.rm=TRUE),
      aum2=aum2011_median[1],
      treated=treated[1],
      accuse_fees=sign(sum(accuse_fees)),accuse_portfolio=sign(sum(accuse_portfolio)),accuse_churning=sign(sum(accuse_churning)),
      accuse_unauth=sign(sum(accuse_unauth)),accuse_frivolous=sign(sum(accuse_frivolous)),accuse_fraid=sign(sum(accuse_fraid)),
      accuse_fiduciary=sign(sum(accuse_fiduciary)),accuse_misrep=sign(sum(accuse_misrep)),accuse_suit=sign(sum(accuse_suit)),accuse_n5000=sign(sum(accuse_n5000)),accuse_n10000=sign(sum(accuse_n10000))
    ),by=list(firmcrd,year,post)]
    firm_panel2[,cl:=state.final]
    
    # wrap analyses as a function:
    r=function(){
      
      felm(data=firm_panel2[year %in% seq(2009,2014),], y(accuse_fiduciary*100/N) ~ 
             post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% 
        summary2
      felm(data=firm_panel2[year %in% seq(2009,2014),], y(accuse_churning*100/N) ~ 
             post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% 
        summary2
      felm(data=firm_panel2[year %in% seq(2009,2014),], y(accuse_suit*100/N) ~ 
             post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% 
        summary2
      felm(data=firm_panel2[year %in% seq(2009,2014),], y(accuse_misrep*100/N) ~ 
             post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% 
        summary2
      
      felm(data=firm_panel2[year %in% seq(2009,2014),], y(accuse_fraid*100/N) ~ 
             post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% 
        summary2
      
      felm(data=firm_panel2[year %in% seq(2009,2014),], y(accuse_n5000*100/N) ~ 
             post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% 
        summary2
      
      felm(data=firm_panel2[year %in% seq(2009,2014),], y(accuse_n10000*100/N) ~ 
             post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% 
        summary2
      
      
      felm(data=firm_panel2[year %in% seq(2009,2014),], y(accuse_fees*100/N) ~ 
             post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% 
        summary2
      
      
      
    }
    y=function(x) sign(x)*100
    l=list(); r() ;         stargazer(l,type='text',omit.stat='ser',no.space=TRUE)
  }
  
}

# Table 5: Regulatory Action----
if(runRegAction){
  
  g=individual_panel
  g[,firmpost:=paste(firmcrd,post)]
  
  l=list()
  felm(data=g %>% subset(year %in% seq(2009,2014)),sign2(complaint.reg) ~ treated*post  |firmcrd+year|0|state.final) %>% summary2
  felm(data=g %>% subset(year %in% seq(2009,2014)),sign2(complaint.reg) ~ treated*post  |firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  felm(data=g %>% subset(year %in% seq(2009,2014)),sign2(complaint.reg) ~ treated*post |year+firmcrd+I(paste(branchclean,post))|0|state.final) %>% summary2
  felm(data=g %>% subset(year %in% seq(2009,2014)),sign2(complaint.reg) ~ treated*post*sign(complaint.lag+complaint.sign) |
         firmcrd+I(paste(branchclean,post))+I(paste(state.final,year))|0|state.final) %>% summary2
  felm(data=g %>% subset(year %in% seq(2009,2014)),sign2(complaint.reg) ~ treated*post*sign(complaint.lag+complaint.sign) |
         firmpost+year|0|state.final) %>% summary2
  a=list(c('firm fe?',grepfe(l,'firm')),
         c('year fe?',grepfe(l,'year')),
         c('branchpost fe?',grepfe(l,'branchclean')),
         c('stateyear fe?',grepfe(l,'state')),
         c('firmpost fe?',grepfe(l,'firmpost'))
  )
  stargazer(l,type='text',omit.stat='ser',no.space=TRUE,
            add.lines=a) 
  
  if (runRegAction_extra) {
    l=list()
    felm(data=g %>% subset(year %in% seq(2009,2014)),sign2(complaint.reg) ~ treated*post*complaint.sign |firmcrd+year|0|state.final) %>% summary2
    felm(data=g %>% subset(year %in% seq(2009,2014)),sign2(complaint.reg) ~ treated*post*complaint.lag |firmcrd+year|0|state.final) %>% summary2
    felm(data=g %>% subset(year %in% seq(2009,2014)),sign2(complaint.reg) ~ treated*post*sign(complaint.lag+complaint.sign) |
           firmcrd+I(paste(branchclean,post))+I(paste(state.final,year))|0|state.final) %>% summary2
    felm(data=g %>% subset(year %in% seq(2008,2015)),sign2(complaint.reg) ~ treated*post*sign(complaint.lag+complaint.sign) |
           firmcrd+I(paste(branchclean,post))+I(paste(state.final,year))|0|state.final) %>% summary2
    felm(data=g %>% subset(year %in% seq(2009,2014)),sign2(complaint.reg) ~ treated*post*sign(complaint.lag+complaint.sign) |
           firmpost+year|0|state.final) %>% summary2
    
    
    felm(data=g %>% subset(year %in% seq(2008,2015)),sign2(complaint.reg) ~ treated*post*sign(complaint.lag+complaint.sign) |
           firmpost+year|0|state.final) %>% summary2
    
    a=list(c('firm fe?',grepfe(l,'firm')),
           c('year fe?',grepfe(l,'year')),
           c('branchpost fe?',grepfe(l,'branchclean')),
           c('stateyear fe?',grepfe(l,'state')),
           c('firmpost?',grepfe(l,'firmpost'))
    )
    stargazer(l,type='text',title='Regulatory action times complaint',omit.stat='ser',no.space=TRUE,
              add.lines=a) ;
  }
  rm(g) ; gc()  
}

# Table 6: Regulator Budget----
if(runBudgetSalary){
    
  x2=merge(individual_panel,budget_1999[,list(staff_per_firm,staff_per_firm0=coalesce(staff_per_firm,0)
                                              ,state.final)],
           by='state.final')
  g=merge(individual_panel,salaries,by=c('state.final'))
  g3=merge(individual_panel,budget_new,by=c('state.final'))
  
  l=list()
  years = seq(2009,2014)
  felm(data=x2 %>% subset(year %in% years),sign2(complaint.sign) ~ treated*post*staff_per_firm |firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  felm(data=x2 %>% subset(year %in% years),sign2(complaint.sign) ~ treated*post*staff_per_firm0 |firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  felm(data=g %>% subset(year %in% years),sign2(complaint.sign) ~ treated*post*log(average_salary) |
         firmcrd+I(paste(state.final,year))| 0 | state.final) %>% summary2
  felm(data=g %>% subset(year %in% years ),sign2(complaint.sign) ~ treated*post*cut(average_salary,breaks = c(0,50000,70000,200000)) |
         firmcrd+I(paste(state.final,year))| 0 | state.final) %>% summary2
  felm(data=g3[year %in% years,],sign2(complaint.sign) ~ post*treated*winsor(budget2012/budget2009,l=0.05,h=0.95)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  felm(data=g3[year %in% years,],sign2(complaint.sign) ~ post*treated*I(budget2012>budget2009)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  felm(data=g3[year %in% years,],sign2(complaint.sign) ~ post*treated*sign(rankz_budget>=0.75)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  
  a=list(c('state fe',grepfe(l,'state')),
         c('firm fe',grepfe(l,'firm')),
         c('year fe',grepfe(l,'year'))
  )
  
  stargazer(l,type='text',omit.stat='ser',no.space = TRUE,
            add.lines=a,keep=c('treated','lag'))
  
  rm(list = c("g3", "g")) ; gc()
  
  # additional results for robustness.
  if (runBudgetSalary_extra) {
    
    l=list()
    years = seq(2009,2014)
    felm(data=x2 %>% subset(year %in% years),sign2(complaint.sign) ~ treated*post*staff_per_firm |firmcrd+year|0|state.final) %>% summary2
    felm(data=x2 %>% subset(year %in% years),sign2(complaint.sign) ~ treated*post*staff_per_firm |firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    felm(data=x2 %>% subset(year %in% years),sign2(complaint.sign) ~ treated*post*staff_per_firm0 |firmcrd+year|0|state.final) %>% summary2
    felm(data=x2 %>% subset(year %in% years),sign2(complaint.sign) ~ treated*post*staff_per_firm0 |firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    a=list(c('state fe',grepfe(l,'state')),
           c('firm fe',grepfe(l,'firm')),
           c('year fe',grepfe(l,'year')),
           c('period',
             ifelse(sapply(l,function(x) x$call %>% as.character %>% unlist %>% grepl(pattern='2009') %>% any),'09/14','08/15'))
    )
    
    stargazer(l,type='text',omit.stat='ser',no.space = TRUE,
              add.lines=a,keep=c('treated','lag'))
    
    rm(x2); gc() # treat your memory
    
    # Salary
    g = merge(individual_panel,salaries,by=c('state.final'))
    options(scipen=10) # so the numbers look less funny
    l=list()
    years=seq(2009,2014)
    felm(data=g %>% subset(year %in% years),sign2(complaint.sign) ~ treated*post*log(average_salary) |
           firmcrd+I(paste(state.final,year))| 0 | state.final) %>% summary2
    felm(data=g %>% subset(year %in% years ),sign2(complaint.sign) ~ treated*post*cut(average_salary,breaks = c(0,60000,120000,200000)) |
           firmcrd+I(paste(state.final,year))| 0 | state.final) %>% summary2
    felm(data=g %>% subset(year %in% years ),sign2(complaint.sign) ~ treated*post*cut(average_salary,breaks = c(0,70000,100000,200000)) |
           firmcrd+I(paste(state.final,year))| 0 | state.final) %>% summary2
    felm(data=g %>% subset(year %in% years ),sign2(complaint.sign) ~ treated*post*cut(average_salary,breaks = c(0,50000,70000,200000)) |
           firmcrd+I(paste(state.final,year))| 0 | state.final) %>% summary2   
    
    stargazer(l,type='text',omit.stat='ser',no.space=TRUE, 
              add.lines=list(c(
                c('firm fe?',grepfe(l,'firmcrd')),
                c('year fe?',grepfe(l,'year')),
                c('state fe?',grepfe(l,'state')
                ))))
    rm(g) ; gc()
    
    # Regulator Budgets
    budget_new[,rankz_budget:=rankz(budget2012/budget2009)]
    
    g3=merge(individual_panel,budget_new,by=c('state.final'))
    l=list()
    years = seq(2009,2014)
    felm(data=g3[year %in% years,],sign2(complaint.sign) ~ post*treated*winsor(budget2012/budget2009,l=0.05,h=0.95)|firmcrd+year|0|state.final) %>% summary2
    felm(data=g3[year %in% years,],sign2(complaint.sign) ~ post*treated*winsor(budget2012/budget2009,l=0.05,h=0.95)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    felm(data=g3[year %in% years,],sign2(complaint.sign) ~ post*treated*I(budget2012>budget2009)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    felm(data=g3[year %in% years,],sign2(complaint.sign) ~ post*treated*rankz_budget|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    felm(data=g3[year %in% years,],sign2(complaint.sign) ~ post*treated*cut(rankz_budget,breaks = c(-0.01,.25,0.5,0.75,1.01)) |firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    felm(data=g3[year %in% years,],sign2(complaint.sign) ~ post*treated*sign(rankz_budget>=0.75)|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    a=list(c('state-year fe',grepfe(l,'state')),
           c('firm fe',grepfe(l,'firm')),
           c('year fe',grepfe(l,'year')),
           c('period',
             ifelse(sapply(l,function(x) x$call %>% as.character %>% unlist %>% grepl(pattern='2008') %>% any),'08/15','09/14'))
    )
    
    stargazer(l,type='text',omit.stat='ser',no.space = TRUE,
              add.lines=a,keep=c('treated','lag'))
    
    rm(l);rm(g3) ; gc()
  }
}

# Table 7: Regulator Distance----
if(runDistance){
  
  g2=merge(individual_panel,distances_individual,by=c('crd','year'))
  g2[!is.na(branchclean),branchpost:=paste(branchclean,post)]
  g2[,cl:=state.final]
  g2[,treated:=sign(treated)][,post:=sign(post)]
  
  g2[,branchpost:=I(paste(branchclean,post))]
  
  years = seq(2009,2014)
  felm(data=g2 %>% subset(year %in% years)  %>% mutate(dist=I(distBranchFromRegulator>=50)),
       sign2(complaint.sign) ~ treated*post*dist
       |firmcrd+year+branchpost|0|cl) %>% summary2
  felm(data=g2 %>% subset(year %in% years) %>% mutate(dist=distRankedWithinState),sign2(complaint.sign) ~ treated*post*dist
       |firmcrd+year+branchpost|0|cl) %>% summary2
  
  felm(data=g2 %>% subset(year %in% years) %>% mutate(dist=log(distBranchFromRegulator+1)),
       sign2(complaint.sign) ~ treated*post*dist+
         treated*post*log(distBranchFromFINRA+1)+treated*post*log(distBranchFromSEC+1)
       |firmcrd+year+branchpost|0|cl) %>% summary2
  
  felm(data=g2 %>% subset(year %in% years) %>% mutate(dist=log(closestRegulator+1)),
       sign2(complaint.sign) ~ treated*post*dist
       |firmcrd+year+branchpost|0|cl) %>% summary2
  
  # for internet distance
  g3=merge(individual_panel,internet,by='branchclean')
  g3=g3[!is.na(fips) & year %in% years,]
  diff=array()[0]
  
  g3[,diff:=total_residential_prov] ; diff=c(diff,'total_residential_prov')    
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
  
  g3[,diff:=total_residential_prov_nbp] ; diff=c(diff,'total_residential_prov_nbp')    
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
  
  a=list(c('firm,year fe?',grepfe(l,'firm')),
         c('branchpost fe?',grepfe(l,'branchpost')))
  stargazer(l,type='text',no.space=TRUE,omit.stat='ser',add.lines=a)      
  
  rm(list = c("g3", "a", "l")) ; gc()
  
  if (runDistance_extra) {  
    # Raw distance
    l=list()
    
    years = seq(2009,2014)
    felm(data=g2 %>% subset(year %in% years) %>% mutate(dist=log(distBranchFromRegulator+1)),
         sign2(complaint.sign) ~ treated*post*dist
         |firmcrd+year+branchpost|0|cl) %>% summary2
    felm(data=g2 %>% subset(year %in% years)  %>% mutate(dist=I(distBranchFromRegulator>=50)),
         sign2(complaint.sign) ~ treated*post*dist
         |firmcrd+year+branchpost|0|cl) %>% summary2
    
    felm(data=g2 %>% subset(year %in% years) %>% mutate(dist=distRankedWithinState),sign2(complaint.sign) ~ treated*post*dist
         |firmcrd+year+branchpost|0|cl) %>% summary2
    
    felm(data=g2 %>% subset(year %in% years),sign2(complaint.sign) ~
           treated*post*cut(distBranchFromRegulator,c(-1,50,250,500,100000))
         |firmcrd+year+branchpost|0|cl) %>% summary2
    
    felm(data=g2 %>% subset(year %in% years) %>% mutate(dist=log(distBranchFromRegulator+1)),
         sign2(complaint.sign) ~ treated*post*dist+
           treated*post*log(distBranchFromFINRA+1)+treated*post*log(distBranchFromSEC+1)
         |firmcrd+year+branchpost|0|cl) %>% summary2
  
    felm(data=g2 %>% subset(year %in% years) %>% mutate(dist=log(closestRegulator+1)),
         sign2(complaint.sign) ~ treated*post*dist
         |firmcrd+year+branchpost|0|cl) %>% summary2
    
    felm(data=g2 %>% subset(year %in% years) %>% mutate(dist=diffStateBranch),sign2(complaint.sign) ~ 
           treated*post*dist
         |firmcrd+year+branchpost|0|cl) %>% summary2
    
  
    a=list(c('firm,year fe?',grepfe(l,'firm')),
           c('branchpost fe?',grepfe(l,'branchpost')))
    stargazer(l,type='text',no.space=TRUE,omit.stat='ser',add.lines=a)      
  
  
  # Changes in distance

  l=list()
  g2[,closerOfFinraSEC:=ifelse(distBranchFromFINRA<=distBranchFromSEC,
                               distBranchFromSEC,distBranchFromSEC)]
  g2[,list(abs(closestRegulator-closerOfFinraSEC) %>% winsor)] %>% summary
  g2[,changeDist:=abs(closestRegulator -      closerOfFinraSEC)]
  g2[,branchpost:=paste(branchclean,post)]
  felm(data=g2 %>% subset(year %in% seq(2009,2014)),sign2(complaint.sign) ~ 
         treated*post*winsor(changeDist,h=0.95,l=0.05)+
         treated*post*log1(closerOfFinraSEC)
       |firmcrd+year+branchpost|0|cl) %>% summary2
  felm(data=g2 %>% subset(year %in% seq(2009,2014)),sign2(complaint.sign) ~ 
         treated*post*cut(changeDist,breaks = c(-1,50,150,300,100000))+
         treated*post*log1(closerOfFinraSEC)
       |firmcrd+year+branchpost|0|cl) %>% summary2
  
  N=50
  g2[,movedCloser50:=(closestRegulator<=50 & distBranchFromSEC>50)]
  g2[,movedCloser250:=(closestRegulator<=250 & distBranchFromSEC>250)]
  felm(data=g2 %>% subset(year %in% seq(2009,2014)),sign2(complaint.sign) ~ 
         treated*post*movedCloser50+ treated*post*log1(closerOfFinraSEC)
       |firmcrd+year+branchpost|0|cl) %>% summary2
  
  felm(data=g2 %>% subset(year %in% seq(2009,2014)),sign2(complaint.sign) ~ 
         treated*post*movedCloser250+ treated*post*log1(closerOfFinraSEC)
       |firmcrd+year+branchpost|0|cl) %>% summary2      
  a=list(c('firm,year fe?',grepfe(l,'firm')),
         c('branchpost fe?',grepfe(l,'branchpost')))
  stargazer(l,type='text',no.space=TRUE,omit.stat='ser',add.lines=a)          
  
  # internet
  g3=merge(individual_panel,internet,by='branchclean')
  g3=g3[!is.na(fips) & year %in% seq(2009,2014),]
  l=list()
  diff=array()[0]
  g3[,diff:=rfc_per_1000_hhs] ; diff=c(diff,'rfc_per_1000_hhs')
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year|0|state.final) %>% summary2
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
  g3[,diff:=I(rfc_per_1000_hhs>=4)] ; diff=c(diff,diff,'I(rfc_per_1000_hhs>=4)')
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
  
  g3[,diff:=total_prov] ; diff=c(diff,'total_prov')    
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
  g3[,diff:=total_residential_prov] ; diff=c(diff,'total_residential_prov')    
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
  g3[,diff:=log(total_residential_prov)] ; diff=c(diff,'log total_residential_prov')    
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
  g3[,diff:=total_residential_prov_nbp] ; diff=c(diff,'total_residential_prov_nbp')    
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
  g3[,diff:=log1(total_residential_prov_nbp)] ; diff=c(diff,'log total_residential_prov_nbp')    
  
  felm(data=g3,sign2(complaint.sign)~post*treated*diff|firmcrd+year+I(paste(fips,post))|0|state.final) %>% 
    summary2
  stargazer(l,type='text',omit.stat='ser',
            add.lines=list(c('firm fe?',rep('Y',length(l))),
                           c('county-post fe?',grepfe(l,'fips')),
                           c('variable',diff %>% gsub('\\_','',.))))
  }
}

# Table 8: Market Power----
  # Note: for individuals
if(runIndivMarketPower){
  # non-accredited investors
  g=merge(client_types,individual_panel,by='firmcrd')
  years = seq(2009,2014)
  g=g[year %in% years,]
  
  # Panel A:
  l=list()
  g %>%
    mutate(individualClients=ifelse(individualClients %in% c('More than 75%','100%','76-99%'),1,0)) %>%
    felm(data=.,sign2(complaint.sign) ~ post*treated*individualClients|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
  g %>%
    mutate(individualClients=ifelse(individualClients %in% c('More than 75%','100%','76-99%'),1,0)) %>%
    felm(data=.,sign2(complaint.sign) ~ post*treated*individualClients|firmcrd+I(paste(state.final,year))+I(paste(branchclean,post))|0|state.final) %>% summary2
  
  g=copy(individual_panel) %>% setDT
  g=merge(g[year==2012,list(numCrds=length(crd),aum2011_median=aum2011_median %>% mean),by=list(firmcrd)] %>% 
            mutate(aumToPeople=winsor(aum2011_median/numCrds)) %>%
            mutate(aumToPeopleRank=rankz(aumToPeople))
          ,g,by='firmcrd')
  setDT(g)
  g[,cl:=state.final]
  felm(data=g[year %in% years,],sign2(complaint.sign) ~ treated*post*log(1+aumToPeople)|firmcrd+I(paste(state.final,year))|0|cl) %>% summary2
  felm(data=g[year %in% years,],sign2(complaint.sign) ~ treated*post*I(aumToPeopleRank>=0.75)|firmcrd+I(paste(year))|0|cl) %>% summary2
  
  stargazer(l,type='text',omit.stat='ser',
            add.lines=list(c('firm fe?',rep('Y',length(l))),
                           c('year fe?', grepfe(l,'year')),
                           c('branch-post fe?',grepfe(l,'branch')),
                           c('state-year fe?',grepfe(l,'state'))))
  
  
  # Panel B:
  # Note: See below for explanation of how HHI is calculated
  x5=merge(branchclean_hhi_employment,individual_panel,by=c('branchclean'))
  x6=individual_panel %>% merge(crosswalk_fips_branch,by='branchclean')
  
  x6=merge(county_hhi_employment,x6,by=c('fips'))
  
  l=list()
  x5=x5[year %in% years & !is.na(branchclean),]
  x6=x6[year %in% years & !is.na(branchclean),]
  
  felm(data=x6,sign2(complaint.sign) ~ post*treated*log(npeople)+
         post*treated*log(hhiPeople*100+1)|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
  felm(data=x5,sign2(complaint.sign) ~ post*treated*log(npeople)+
         post*treated*log(hhiPeople*100+1)|firmcrd+year+I(paste(branchclean,post))|0|state.final) %>% summary2
  
  # AUM-based HHI for State
  merge(hhi,individual_panel,by='state.final') %>% 
    subset(year %in% seq(2009,2014)) %>%
    felm(data=.,sign2(complaint.sign) ~ post*treated*log(top4AUM)|firmcrd+
           I(paste(state.final,year))|0|state.final) %>% summary2
  
  merge(hhi,individual_panel,by='state.final') %>% 
    subset(year %in% seq(2009,2014)) %>%
    felm(data=.,sign2(complaint.sign) ~ post*treated*log(hhiAUM)|firmcrd+I(paste(state.final,year))|0|state.final) %>% 
    summary2
  
  stargazer(l,type='text',omit.stat='ser',
            add.lines=list(c('firm fe?',rep('Y',length(l))),
                           c('year fe?', grepfe(l,'year')),
                           c('branch-post fe?',grepfe(l,'branch')),
                           c('state-year fe?',grepfe(l,'state')),
                           c('aggregation?','county','city',rep('',3))))
  
  
  if (runIndivMarketPower_extra) {
  {
    g=merge(client_types,individual_panel,by='firmcrd')
    g=g[year %in% years,]
    l=list()
    g %>% 
      mutate(individualClients=ifelse(individualClients %in% c('More than 75%','100%','76-99%'),1,0)) %>%
      felm(data=.,sign2(complaint.sign) ~ post*treated*individualClients|firmcrd+year|0|state.final) %>% summary2
    g %>%
      mutate(individualClients=ifelse(individualClients %in% c('More than 75%','100%','76-99%'),1,0)) %>%
      felm(data=.,sign2(complaint.sign) ~ post*treated*individualClients|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary2
    
    g %>%
      mutate(individualClients=ifelse(individualClients %in% c('More than 75%','100%','76-99%'),1,0)) %>%
      felm(data=.,sign2(complaint.sign) ~ post*treated*individualClients|firmcrd+I(paste(state.final,year))|0|state.final) %>% summary
    g %>%
      mutate(individualClients=ifelse(individualClients %in% c('More than 75%','100%','76-99%'),1,0)) %>%
      felm(data=.,sign2(complaint.sign) ~ post*treated*individualClients|firmcrd+I(paste(state.final,year))+I(paste(branchclean,post))|0|state.final) %>% summary2
    
    stargazer(l,type='text',omit.stat='ser',
              add.lines=list(c('firm fe?',rep('Y',length(l))),
                             c('branch-post fe?',grepfe(l,'branch')),
                             c('state-year fe?',grepfe(l,'state'))))
  }
  
  # aum to people rank
  {
    
    l=list()
    rm(g) ; gc()
    g=copy(individual_panel) %>% setDT
    g=merge(g[year==2012,list(numCrds=length(crd),aum2011_median=aum2011_median %>% mean),by=list(firmcrd)] %>% 
              mutate(aumToPeople=winsor(aum2011_median/numCrds)) %>%
              mutate(aumToPeopleRank=rankz(aumToPeople))
            ,g,by='firmcrd')
    setDT(g)
    g[,cl:=state.final]
    l=list()
    felm(data=g[year %in% years,],sign2(complaint.sign) ~ treated*post*aumToPeopleRank|firmcrd+year|0|cl) %>% summary2
    felm(data=g[year %in% years,],sign2(complaint.sign) ~ treated*post*log(1+aumToPeople)|firmcrd+year|0|cl) %>% summary2
    felm(data=g[year %in% years,],sign2(complaint.sign) ~ treated*post*I(aumToPeopleRank>=0.75)|firmcrd+year|0|cl) %>% summary2
    
    felm(data=g[year %in% years,],sign2(complaint.sign) ~ treated*post*log(1+aumToPeople)|firmcrd+I(paste(state.final,year))|0|cl) %>% summary2
    felm(data=g[year %in% years,],sign2(complaint.sign) ~ treated*post*I(aumToPeopleRank>=0.75)|firmcrd+I(paste(state.final,year))|0|cl) %>% summary2
    
    stargazer(l,type='text',omit.stat='ser',add.lines=list(c('firm+year fe',grepfe(l,'firmcrd')),
                                                           c('statebyyear fe',grepfe(l,'state'))
                                                           )) 
  }  
  
  # bargaining power - aum of the state concentration
  {
    
    l=list()
    merge(hhi,individual_panel,by='state.final') %>% 
      subset(year %in% seq(2009,2014)) %>%
      felm(data=.,sign2(complaint.sign) ~ post*treated*log(top4AUM)|firmcrd+
             I(paste(state.final,year))|0|state.final) %>% summary2
    
    merge(hhi,individual_panel,by='state.final') %>% 
      subset(year %in% seq(2009,2014)) %>%
      felm(data=.,sign2(complaint.sign) ~ post*treated*log(hhiAUM)|firmcrd+I(paste(state.final,year))|0|state.final) %>% 
      summary2
    
    stargazer(l,type='text',omit.stat='ser',add.lines=list(c('firm+year fe',grepfe(l,'firmcrd')),
                                                           c('statebyyear fe',grepfe(l,'state'))
    )) 
    
  }
  
  # Employment HHI is based on the full county distribution and the full branch distribution
    # Since the branch is just a city name, it's not a real aggregation unit
    # The county is as well
  {
    
    x5=merge(branchclean_hhi_employment,individual_panel,by=c('branchclean'))
    x6=individual_panel %>% merge(crosswalk_fips_branch,by='branchclean')
    
    x6=merge(county_hhi_employment,x6,by=c('fips'))
    
    l=list()
    x5=x5[year %in% seq(2009,2014) & !is.na(branchclean),]
    x6=x6[year %in% seq(2009,2014) & !is.na(branchclean),]
    felm(data=x5,sign2(complaint.sign) ~ post*treated*log(npeople)+
           post*treated*hhiPeople|firmcrd+year+I(paste(branchclean,post))|0|state.final) %>% summary2
    
    felm(data=x5,sign2(complaint.sign) ~ post*treated*log(npeople)+
           post*treated*log(hhiPeople*100+1)|firmcrd+year+I(paste(branchclean,post))|0|state.final) %>% summary2

    felm(data=x6,sign2(complaint.sign) ~ 
           post*treated*hhiPeople|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2    
    felm(data=x6,sign2(complaint.sign) ~ post*treated*log(npeople)+
           post*treated*hhiPeople|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
    
    felm(data=x6,sign2(complaint.sign) ~ post*treated*log(npeople)+
           post*treated*log(hhiPeople*100+1)|firmcrd+year+I(paste(fips,post))|0|state.final) %>% summary2
    
    stargazer(l,type='text',omit.stat='ser',no.space=TRUE,
              add.lines=list(c('firm fe?',rep('Y',length(l))),
                             c('county-post fe?',length(l)),
                             c('aggregation?','city','city','county','county','county')
              ),keep=':')
    
    
  }
  
  }
}

# Table 9: Recidivism----
if(runRepeatComplaints){
  history=read_fst('data_for_distribution/complaints_before_2010.fst') %>% setDT
  g=merge(individual_panel,history,by='crd')
  l=list()
  felm(data=g[ year %in% seq(2009,2014) & aum2011_median<=(10^9),],I(100*complaint.sign) ~ 
         sign(history_pre2010)*treated*post
       |I(paste(firmcrd)) +year|0|state.final) %>% summary2
  felm(data=g[ year %in% seq(2009,2014) & aum2011_median<=(10^9),],I(100*complaint.sign) ~ 
         sign(history_pre2010)*treated*post
       |I(paste(firmcrd,post)) +year|0|state.final) %>% summary2
  
  felm(data=g[ year %in% seq(2009,2014) & aum2011_median<=(3*10^8),],I(100*complaint.sign) ~ 
         sign(history_pre2010)*treated*post+treated*year
       |I(paste(firmcrd,post)) +year|0|state.final) %>% summary2
  
  stargazer(l,type='text',omit.stat='ser',no.space = TRUE,
            add.lines=list(c('fe?','firm + year','firm*post+year','firm*post+year'),
                           c('sample','full','full','<=300mm')),
            notes = c("Column (3) not shown in the main paper but included here as an additional robustness check."));
}

# Table 10: Dual Registered----
if(runDualRegistered){
  {
    # dual registered analysis
    {
      
      
      firm_panel2=firm_panel %>% copy %>% merge(fract_dual_registered_firm,by=c('firmcrd','year'))
      firm_panel2[year %in% seq(2009,2014),]$fractEverDualRegister %>% summary
      # the median is 0.2857
      firm_panel2[,post:=sign(year>=2012)]
      
      l=list()
      felm(data=firm_panel2[fractEverDualRegister>=0.2857,]  %>% subset(year %in% seq(2009,2014)),
           sign2(n_complaints)~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% summary2
      felm(data=firm_panel2[fractEverDualRegister<=0.2857,]  %>% subset(year %in% seq(2009,2014)),
           sign2(n_complaints)~ post*treated+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% summary2
      
      felm(data=firm_panel2[,] %>% subset(year %in% seq(2009,2014)),
           sign2(n_complaints)~ post*treated*I(fractEverDualRegister>=0.2857)+log(N)+I(log(N)^2)+I(log(N)^3)|
             firmcrd+year|0|state.final) %>% summary2
      x5=inner_join(individual_panel,dual_registered_indiv) 
      felm(data=x5[,] %>% subset(year %in% seq(2009,2014)),
           sign2(complaint.sign)~ post*treated*dualRegister |
             firmcrd+year|0|state.final) %>% summary2
      
      
      stargazer(l,type='text',no.space=TRUE,omit.stat='ser')
      
    }
    
    # broker dealer
    if (runDualRegistered_extra)
    {
      
      firm_panel2=copy(firm_panel)
      firm_panel2[,brokerDealer:=firmcrd %in% bdCRDs]
      sign100=sign2
      
      l=list()
      felm(data=firm_panel2[year %in% seq(2009,2014),],
           sign100(n_complaints)~ post*treated*brokerDealer+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% summary2
      
      felm(data=firm_panel2[year %in% seq(2009,2014) & brokerDealer==TRUE,],
           sign100(n_complaints)~ post*treated*brokerDealer+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% summary2
      
      felm(data=firm_panel2[year %in% seq(2009,2014) & brokerDealer==FALSE,],
           sign100(n_complaints)~ post*treated*brokerDealer+log(N)+I(log(N)^2)+I(log(N)^3)|firmcrd+year|0|state.final) %>% summary2
      
      stargazer(l,type='text',omit.stat='ser',no.space=TRUE)
      
    }
    
  }
}

# Figures: Parallel Trend Plots----
if(runPTgraphs){
  
  controlLabel='SEC registered '
  
  routine=function(out1,label=controlLabel,numse=1){
    require(broom)      
    toPlot1 = tidy(out1) %>% mutate(year = removeNotNumbers(term) %>% as.numeric,
                                    treated = !grepl(term, pattern = "FALSE")) %>% 
      filter(!is.na(year)) %>% mutate(variable = ifelse(treated, yes = "Switch to State Registered  ", 
                                                        no = controlLabel),
                                      value = estimate) %>%
      arrange(year) %>% group_by(treated) %>%
      mutate(value = value - mean(value, na.rm = T));
    toPlot1 %>% setDT
    toPlot1=toPlot1[year>=2005,]
    
    plotPT = ggplot(toPlot1) +
      geom_line(aes(x = year, y = value, col = variable, linetype = variable), size = 1) +
      scale_shape_manual(values=c(16,2))+ geom_point(aes(x = year, y = value, shape = variable, col = variable), size = 5) + 
      geom_errorbar(aes(x = year, ymin = value-std.error*numse, ymax = value+std.error*numse, col = variable), width = 0.15) +
      geom_vline(xintercept = 2012+1/12, linetype = "dashed", alpha = 0.5, size = 1.25, col = "blue") + 
      labs(title = "Complaints across RIAs", y = "P(Complaints)", x = "", subtitle = "") +
      scale_x_continuous(breaks = seq(min(toPlot1$year),max(toPlot1$year)),
                         labels = seq(min(toPlot1$year),max(toPlot1$year)))    +
      theme_BC() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),         
                         panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
      scale_color_discrete() + theme(legend.key.width = unit(2,"cm"));
    
    print(plotPT)
  }
  
  firm_panel[ year %in% seq(2006,2015),] %>% felm(data=.,sign2(n_complaints) ~ 
                                                    as.factor(year):I(treated==TRUE)+log(N)+I(log(N)^2)+I(log(N)^3)-1|0|
                                                    0|state.final)  -> out1
  routine(out1,numse=1)
  
  
  firm_panel[ year %in% seq(2006,2015) & aum2<=(10^8*3),] %>% felm(data=.,sign2(n_complaints) ~ 
                                                                     as.factor(year):I(treated==TRUE)+log(N)+I(log(N)^2)+I(log(N)^3)-1|0|
                                                                     0|state.final)  -> out1
  routine(out1,numse=1)
  
  
  firm_panel[ year %in% seq(2006,2015)  & aum2<=(10^8),] %>% felm(data=.,sign2(n_complaints) ~ 
                                                                    as.factor(year):I(treated==TRUE)+log(N)+I(log(N)^2)+I(log(N)^3)-1|0|
                                                                    0|state.final)  -> out1
  routine(out1,numse=2)
  
  # Plot 4
  output=function(haha,n=1.645){
    hi=haha  %>% broom::tidy() %>% setDT
    hi[grepl(term,pattern=':treated') | grepl(term,pattern='treated:') | grepl(term,pattern='treatedTRUE:'),] %>% mutate(year=removeNotNumbers(`term`) %>% as.integer) %>% 
      mutate(se=`std.error`) %>%
      ggplot(aes(x=year,y=estimate)) +
      geom_point(size=3) +
      geom_line(linetype='dashed') +
      geom_hline(yintercept = 0, linetype = "dashed", col = "blue", alpha = 0.5) +
      geom_vline(xintercept = 2012, linetype = "dashed", col = "black", alpha = 0.5) +
      geom_errorbar(aes(ymin=estimate-n*se,ymax=estimate+n*se,alpha=0.5),color='black',width = 0.15) +
      ggthemes::theme_stata() +
      theme(legend.position = "none") %>% print
  }
  
  
  felm(data=firm_panel[year %in% seq(2008,2015) & aum2<=(1000*10^6),],sign2(n_complaints) ~ treated*as.factor(year)+log(N)+
         I(log(N)^2)+I(log(N)^3)|firmcrd+I(paste(state.final,year))|0|state.final) %>% output(1.65) +
    theme_bw()
  
  
  
}


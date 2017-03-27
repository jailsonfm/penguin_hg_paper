
#*** Manuscript in preparation - Mercury in Penguins

# ------ Libraries and Functions --------------------------------- -----
library(reshape)
library(agricolae)
library(stats)
library(multcomp)
library(pastecs)
library(car)
library(MASS) 
library(ggplot2) 
library(ellipse) 
library(plyr) # detach("package:plyr", unload=TRUE)
library(FactoMineR) # detach("package:Rmisc", unload=TRUE)
library(mgcv)
library(gridExtra)  #  multiplot ggplot2
library(dplyr)  #  multiplot ggplot2 
tukey <- function(data, group, method=c("Tukey", "Games-Howell")){
        OK <- complete.cases(data, group)                       
        data <- data[OK]
        group <- factor(group[OK])
        n <- tapply(data, group, length)                        
        a <- length(n)                                          
        phi.e <- sum(n)-a                                       
        Mean <- tapply(data, group, mean)                        
        Variance <- tapply(data, group, var)                     
        result1 <- cbind(n, Mean, Variance)                      
        rownames(result1) <- paste("Group", 1:a, sep="")
        method <- match.arg(method)
        if (method == "Tukey") {                                 
                v.e <- sum((n-1)*Variance)/phi.e                 
                t <- combn(a, 2, function(ij)                    
                                        abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
                p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)       
                Tukey <- cbind(t, p)                                     
                rownames(Tukey) <- combn(a, 2, paste, collapse=":")
                return(list(result1=result1, Tukey=Tukey, phi=phi.e, v=v.e))
        }
        else {                                                   
                t.df <- combn(a, 2, function(ij) {               
                                        t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
                                        df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
                                        return(c(t, df))} )
                t <- t.df[1,]
                df <- t.df[2,]
                p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)  
                Games.Howell <- cbind(t, df, p)                  
                rownames(Games.Howell) <- combn(a, 2, paste, collapse=":")
                return(list(result1=result1, Games.Howell=Games.Howell))
        }
} # Games-Howell test
#Source: http://www.psych.yorku.ca/cribbie/6130/games_howell.R


# ------ Data Entry and organization ----------------------------- -----

setwd("/Users/jfu/Documents/Projects & Manuscripts/manuscripts/Mercury Penguins/Data_&_Rcode")
dir()
peng.dat <- read.csv("/Users/jfu/Documents/Projects & Manuscripts/manuscripts/Mercury Penguins/Data_&_Rcode/data_peng_ok.csv", sep=",", header=TRUE)
#peng.dat <- rename(penguin, c(Length = "Length", Year="Year", Code="Code", Weight="Weight"))
str(peng.dat)   ;  names(peng.dat)
#peng[peng.dat == 0.000] <- NA
peng.dat$Year <- factor(peng.dat$Year)
head(peng.dat)  ;  tail(peng.dat)
#

peng_year_rj <- filter(peng.dat, Area == "Rio") # Filtering the data for comparison among year in the Rio data
peng_area_2008 <- filter(peng.dat, Year=="2008") # Selecting data just for 2008 for the comparizon Areas

resu <- peng_area_2008 %>%
        group_by(Area) %>%
        summarise(n = n(), Mean = mean(Hg_L), se = sd(Hg_L)/sqrt(n), sd = sd(Hg_L))

peng_year_rj$body_cond <- peng_year_rj$Weight/peng_year_rj$Length

# ------ BoxPlot Hg and biometry across Years for Rio ------------ -----
year.HgM <- ggplot(data = peng_year_rj, aes(x = Year, y = Hg_M)) + 									
        geom_boxplot(aes(fill=factor(Year)),fill="#66A1D2", colour="#27258D", alpha=0.45) +			
        theme(panel.background = element_rect(fill='white', colour='black'), 
              axis.text.x = element_blank(), axis.title.x = element_blank()) +
        ylab("Hg Level (in Muscle)")+ xlab("Years") 	+	
        theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  											# + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)

year.HgL <- ggplot(data = peng_year_rj, aes(x = Year, y = Hg_L)) + 	
        geom_boxplot(aes(fill=factor(Year)),fill="#66A1D2", colour="#27258D", alpha=0.45)+	
        theme(panel.background = element_rect(fill='white', colour='black'), 
              axis.text.x = element_blank(), axis.title.x = element_blank()) +
        ylab("Hg Level (in Liver)")+
        theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  											# + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)

Year.Len <- ggplot(data = peng_year_rj, aes(x = Year, y = Length)) + 	
        geom_boxplot(aes(fill=factor(Year)),fill="#66A1D2", colour="#27258D", alpha=0.45)+	
        theme(panel.background = element_rect(fill='white', colour='black'), 
              axis.text.x  = element_text(vjust=0.5, size=12, color="black"), axis.title.x  = element_blank()) +	
        ylab("Length (in cm)")	+			
        theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),legend.position="none")  											# + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)

Year.Weig <- ggplot(data = peng_year_rj, aes(x = Year, y = Weight)) + 		
        geom_boxplot(aes(fill=factor(Year)),fill="#66A1D2", colour="#27258D", alpha=0.45)+			
        theme(panel.background = element_rect(fill='white', colour='black'), 
              axis.text.x  = element_text(vjust=0.5, size=12, color="black"), axis.title.x  = element_blank())  +
        ylab("Weight (in Kg)")+ xlab("Years") 	+	
        theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  									# + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)

Year.body_cond <- ggplot(data = peng_year_rj, aes(x = Year, y = body_cond)) + 		
        geom_boxplot(aes(fill=factor(Year)),fill="#66A1D2", colour="#27258D", alpha=0.45)+			
        theme(panel.background = element_rect(fill='white', colour='black'), 
              axis.text.x  = element_text(vjust=0.5, size=12, color="black"), axis.title.x  = element_blank())  +
        ylab("body_cond (in Kg)")+ xlab("Years") 	+	
        theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none") 


peng_year_rj$body_cond
# png(filename="/Users/jfu/Documents/Projects & Manuscripts/manuscripts/Mercury Penguins/boxplot_ggplot.png", width=8.5, height=7, res=1500, units="in")
 grid.arrange(year.HgM, year.HgL, Year.Len, Year.Weig, ncol = 2)
#dev.off()

# ------ ANOVA for Biometry and Hg across the years for Rio  --- -----

#-------------------- Welch's F test Hg Muscle across Years ------
 
 
by(peng_year_rj$Hg_M , peng_year_rj$Year, stat.desc)   #  library("pastecs")
leveneTest(peng_year_rj$Hg_M , peng_year_rj$Year, center=mean)  #  library(car)  Signif! Unequal variance!
bartlett.test(Hg_M ~ Year, data= peng_year_rj)  # Signif! Unequal variance!
muscHg.year <- oneway.test(Hg_M ~ Year, data = peng_year_rj); muscHg.year   #  library(stats)  ; Welch's F Test.
tukey(peng_year_rj$Hg_M , peng_year_rj$Year, method="Games-Howell") # Games-Howell Post Hoc Test in R

#-------------------- Welch's F test Hg Liver across Years ---
by(peng_year_rj$Hg_L , peng_year_rj$Year, stat.desc)   #  library("pastecs")
leveneTest(peng_year_rj$Hg_L , peng_year_rj$Year, center=mean)  #  library(car)  Signif! Unequal variance!
bartlett.test(Hg_M ~ Year, data= peng_year_rj)  # Signif! Unequal variance!
livHg.year <- oneway.test(Hg_L ~ Year, data = peng_year_rj) ; livHg.year  #  library(stats)  ; Welch's F Test .
tukey(peng_year_rj$Hg_L , peng_year_rj$Year, method="Games-Howell") # Games-Howell Post Hoc Test in R

#-------------------- Welch's F test Weight values  across Years ---
by(peng_year_rj$Weight , peng_year_rj$Year, stat.desc)   #  library("pastecs")
leveneTest(peng_year_rj$Weight , peng_year_rj$Year, center=median)  #  library(car)  Signif! Unequal variance!
bartlett.test(Weight ~ Year, data= peng_year_rj)  # Signif! Unequal variance!
weight.year <- oneway.test(Weight ~ Year, data = peng_year_rj); weight.year   #  library(stats)  ; Welch's F Test .
tukey(peng_year_rj$Weight , peng_year_rj$Year, method="Games-Howell") # Games-Howell Post Hoc Test in R

#-------------------- Welch's F test Length values  across Years ---
by(peng_year_rj$Length , peng_year_rj$Year, stat.desc)   #  library("pastecs")
leveneTest(peng_year_rj$Length , peng_year_rj$Year, center=mean)  #  library(car)  Signif! Unequal variance!
bartlett.test(Length ~ Year, data= peng_year_rj)  # Signif! Unequal variance!
weight.year <- oneway.test(Length ~ Year, data = peng_year_rj); weight.year   #  library(stats)  ; Welch's F Test .
tukey(peng_year_rj$Length , peng_year_rj$Year, method="Games-Howell") # Games-Howell Post Hoc Test in R

#-------------------- Welch's F test Length values  across Years ---
by(peng_year_rj$body_cond , peng_year_rj$Year, stat.desc)   #  library("pastecs")
leveneTest(peng_year_rj$body_cond , peng_year_rj$Year, center=mean)  #  library(car)  Signif! Unequal variance!
bartlett.test(body_cond ~ Year, data= peng_year_rj)  # Signif! Unequal variance!
body_cond.year <- oneway.test(body_cond ~ Year, data = peng_year_rj); body_cond.year   #  library(stats)  ; Welch's F Test .
tukey(peng_year_rj$body_cond , peng_year_rj$Year, method="Games-Howell") # Games-Howell Post Hoc Test in R

#___________________________________________________________________________________________________________________________________________________________________________

                        #-------------------- ANOVA F test Hg Muscle across Years ---
 			aov_HgM <- aov(Hg_M ~ Year, data = peng_year_rj)  ;    summary(aov_HgM) ##
			outHgM  <- glht(aov_HgM, linfct=mcp(Year="Tukey")); summary(outHgM)   #   library(multcomp) 
			#-------------------- ANOVA F test Hg Liver across Years ---
			 aov_HgL <- aov(Hg_L ~ Year, data = peng_year_rj)  ;    summary(aov_HgL)  ## Anova test Hg_Liver by Years
 			outHgL  <- glht(aov_HgL, linfct=mcp(Year="Tukey")); summary(outHgL)    #   library(multcomp)
 			#-------------------- ANOVA F test Weight across Years ---
 			aov_W <- aov(Weight ~ Year, data = peng_year_rj)  ;    summary(aov_W)
			outW  <- glht(aov_W, linfct=mcp(Year="Tukey")); summary(outHgL)    #   library(multcomp) 
			#-------------------- NOVA F test Length across Years ---
 			aov_TL <- aov(Length ~ Year, data = peng_year_rj)  ;    summary(aov_W)
			outTL  <- glht(aov_TL, linfct=mcp(Year="Tukey")); summary(outTL)   #   library(multcomp) 

# ------ BoxPlot Hg and biometry across Years for Areas ---------- -----
			
Area.HgM <- ggplot(data = peng_area_2008, aes(x = Area, y = Hg_M)) + 									
        geom_boxplot(aes(fill=factor(Area)),fill="#66A1D2", colour="#27258D", alpha=0.45) +			
        theme(panel.background = element_rect(fill='white', colour='black'), 
                axis.text.x = element_blank(), axis.title.x = element_blank()) +
        ylab("Hg Level (in Muscle)")+ xlab("Areas") 	+	
        theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  											# + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)
			
Area.HgL <- ggplot(data = peng_area_2008, aes(x = Area, y = Hg_L)) + 	
        geom_boxplot(aes(fill=factor(Area)),fill="#66A1D2", colour="#27258D", alpha=0.45)+	
        theme(panel.background = element_rect(fill='white', colour='black'), 
                axis.text.x = element_blank(), axis.title.x = element_blank()) +
        ylab("Hg Level (in Liver)")+
        theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  											# + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)
			
Area.Len <- ggplot(data = peng_area_2008, aes(x = Area, y = Length)) + 	
        geom_boxplot(aes(fill=factor(Area)),fill="#66A1D2", colour="#27258D", alpha=0.45)+	
        theme(panel.background = element_rect(fill='white', colour='black'), 
                axis.text.x  = element_text(vjust=0.5, size=12, color="black"), axis.title.x  = element_blank()) +	
        ylab("Length (in cm)")	+			
        theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),legend.position="none")  											# + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)
			
Area.Weig <- ggplot(data = peng_area_2008, aes(x = Area, y = Weight)) + 		
        geom_boxplot(aes(fill=factor(Area)),fill="#66A1D2", colour="#27258D", alpha=0.45)+			
        theme(panel.background = element_rect(fill='white', colour='black'), 
                axis.text.x  = element_text(vjust=0.5, size=12, color="black"), axis.title.x  = element_blank())  +
        ylab("Weight (in Kg)")+ xlab("Areas") 	+	
        theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  									# + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)
			
# png(filename="/Users/jfu/Documents/Projects & Manuscripts/manuscripts/Mercury Penguins/boxplot_ggplot.png", width=8.5, height=7, res=1500, units="in")
grid.arrange(Area.HgM, Area.HgL, Area.Len, Area.Weig, ncol = 2)
#dev.off()

			

# ------ ANOVA for Biometry and Hg across the Areas ------------ -----

#-------------------- Welch's F test Hg Muscle across Areas ---
by(peng_area_2008$Hg_M , peng_area_2008$Area, stat.desc)   #  library("pastecs")
leveneTest(peng_area_2008$Hg_M , peng_area_2008$Area, center=mean)  #  library("car")  Signif! Unequal variance!
bartlett.test(Hg_M ~ Area, data= peng_area_2008)  # Signif! Unequal variance!
muscHg.Area <- oneway.test(Hg_M ~ Area, data = peng_area_2008); muscHg.Area   #  library(stats)  ; Welch's F Test.
tukey(peng_area_2008$Hg_M , peng_area_2008$Area, method="Games-Howell") # Games-Howell Post Hoc Test in R

#-------------------- Welch's F test Hg Liver across Areas ---
by(peng_area_2008$Hg_L , peng_area_2008$Area, stat.desc)   #  library("pastecs")
leveneTest(peng_area_2008$Hg_L , peng_area_2008$Area, center=mean)  #  library(car)  Signif! Unequal variance!
bartlett.test(Hg_M ~ Area, data= peng_area_2008)  # Signif! Unequal variance!
livHg.Area <- oneway.test(Hg_L ~ Area, data = peng_area_2008) ; livHg.Area  #  library(stats)  ; Welch's F Test .
tukey(peng_area_2008$Hg_L , peng_area_2008$Area, method="Games-Howell") # Games-Howell Post Hoc Test in R

#-------------------- Welch's F test Weight values  across Areas ---
by(peng_area_2008$Weight , peng_area_2008$Area, stat.desc)   #  library("pastecs")
leveneTest(peng_area_2008$Weight , peng_area_2008$Area, center=median)  #  library(car)  Signif! Unequal variance!
bartlett.test(Weight ~ Area, data= peng_area_2008)  # Signif! Unequal variance!
weight.Area <- oneway.test(Weight ~ Area, data = peng_area_2008); weight.Area   #  library(stats)  ; Welch's F Test .
tukey(peng_area_2008$Weight , peng_area_2008$Area, method="Games-Howell") # Games-Howell Post Hoc Test in R

#-------------------- Welch's F test Length values  across Areas ---
by(peng_area_2008$Length , peng_area_2008$Area, stat.desc)   #  library("pastecs")
leveneTest(peng_area_2008$Length , peng_area_2008$Area, center=mean)  #  library(car)  Signif! Unequal variance!
bartlett.test(Length ~ Area, data= peng_area_2008)  # Signif! Unequal variance!
weight.Area <- oneway.test(Length ~ Area, data = peng_area_2008); weight.Area   #  library(stats)  ; Welch's F Test .
tukey(peng_area_2008$Length , peng_area_2008$Area, method="Games-Howell") # Games-Howell Post Hoc Test in R
#___________________________________________________________________________________________________________________________________________________________________________

#-------------------- ANOVA F test Hg Muscle across Areas ---
aov_HgM <- aov(Hg_M ~ Area, data = peng_area_2008)  ;    summary(aov_HgM) ##
outHgM  <- glht(aov_HgM, linfct=mcp(Area="Tukey")); summary(outHgM)   #   library(multcomp) 
#-------------------- ANOVA F test Hg Liver across Areas ---
aov_HgL <- aov(Hg_L ~ Area, data = peng_area_2008)  ;    summary(aov_HgL)  ## Anova test Hg_Liver by Areas
outHgL  <- glht(aov_HgL, linfct=mcp(Area="Tukey")); summary(outHgL)    #   library(multcomp)
#-------------------- ANOVA F test Weight across Areas ---
aov_W <- aov(Weight ~ Area, data = peng_area_2008)  ;    summary(aov_W)
outW  <- glht(aov_W, linfct=mcp(Area="Tukey")); summary(outHgL)    #   library(multcomp) 
#-------------------- NOVA F test Length across Areas ---
aov_TL <- aov(Length ~ Area, data = peng_area_2008)  ;    summary(aov_W)
outTL  <- glht(aov_TL, linfct=mcp(Area="Tukey")); summary(outTL)   #   library(multcomp) 

# ------ GLM Hg vs Length & Weight All Years --------------------- ------------
# Check Normality #
qqnorm(peng.dat$Length)      #   hist(peng.dat$Length)
qqnorm(peng.dat$Weight)     #   shapiro.test(peng.dat$Weight)
qqnorm(peng.dat$Hg_L)         #   hist(peng.dat$Hg_L)
qqnorm(peng.dat$Hg_M)       #   hist(peng.dat$Hg_M)
# log normal, exponention, gamma*, 

#---- TRANSFOIRMED DATA log  (Hg Muscle & Liver )
peng.dat$log.HgM <- log(peng.dat$Hg_M)
qqnorm(peng.dat$log.HgM)
peng.dat$log.HgL <- log(peng.dat$Hg_L)
qqnorm(peng.dat$log.HgL)
    
glm.HgM.Len <- glm(Hg_M ~ Length, data = peng_area_2008, family=Gamma(link="log"))   #  GLM Hg Muscle vs Length
		summary(glm.HgM.Len)   #   plot(glm.HgM.Len)
        #  glm.HgM.Len <- glm(Hg_M ~ Length, data = peng.dat, family=gaussian(link="log"))
        #  glm.HgM.Len <- glm(Hg_M ~ Length, data = peng.dat, family=Gamma(link="log"))

glm.HgM.Weig <- glm(Hg_M ~ Weight, data = peng_area_2008, family=Gamma(link="log"))   #  GLM Hg Muscle vs Weight
		summary(glm.HgM.Weig)		#  plot(glm.HgM.Weig)
        #  glm.HgM.Weig <- glm(Hg_M ~ Weight, data = peng.dat, family=gaussian(link="log")) 
        #  glm.HgM.Weig <- glm(Hg_M ~ Weight, data = peng.dat, family=Gamma(link="log")) 
		
glm.HgL.Len <- glm(Hg_L ~ Length, data = peng_area_2008, family=Gamma(link="log"))   #  GLM Hg Liver vs Length
		summary(glm.HgM.Len)  #  plot(glm.HgL.Len)
        #   glm.HgL.Len <- glm(Hg_L ~ Length, data = peng.dat, family=gaussian(link="log"))
        #   glm.HgL.Len <- glm(Hg_L ~ Length, data = peng.dat, family=Gamma(link="log"))  

glm.HgL.Weig <- glm(Hg_L ~ Weight, data = peng_area_2008, family=Gamma(link="log"))   #  GLM Hg Liver vs Weight
		summary(glm.HgL.Weig)#      plot(glm.HgL.Weig)
        #   glm.HgL.Weig <- glm(Hg_L ~ Weight, data = peng.dat, family=gaussian(link="log"))
        #   glm.HgL.Weig <- glm(Hg_L ~ Weight, data = peng.dat, family=Gamma(link="log"))
str(peng.dat)

peng_area_2008

# ------ Plot GLM with ggplot2 ----------------------------------- ----
peng.dat$log.hgm <- peng.dat$Hg_M

# ----- plot Hg_Muscle vs Weight
Hg_Mus_weight <- ggplot(peng.dat, aes(x=Length, y=Hg_L)) + 
        stat_smooth(method = "glm",fill = "green", size = 2)+
      geom_point(shape=21, fill="red", size=2, color="#215A89") +
        theme(panel.border = element_rect(fill=NA,color="black", size=2.0,linetype="solid"),legend.position="none") + 				
        theme(panel.background = element_rect(fill='white', colour='gray')) +														
        ylab("Hg Level (in Muscle)") + xlab("Weight (kg)")   

# ----- plot Hg_Muscle vs Length
Hg_Mus_length <- ggplot(peng.dat, aes(x=Length, y=log(Hg_M))) + 
        stat_smooth(method = "glm", fill = "green", size = 2)+
        geom_point(shape=21, fill="red", size=2, color="#215A89") +
        theme(panel.border = element_rect(fill=NA,color="black", size=2.0,linetype="solid"),legend.position="none") + 				
        theme(panel.background = element_rect(fill='white', colour='gray')) +														
        ylab("Hg Level (in Muscle)") + xlab("Length (cm)")   

# ----- plot Hg_Liver vs Weight
Hg_Liv_length <- ggplot(peng.dat, aes(x=Weight, y=Hg_L)) + 
        stat_smooth(method = "glm", family=Gamma(link="log"), fill = "green", size = 2)+
        geom_point(shape=21, fill="red", size=2, color="#215A89") +
        theme(panel.border = element_rect(fill=NA,color="black", size=2.0,linetype="solid"),legend.position="none") + 				
        theme(panel.background = element_rect(fill='white', colour='gray')) +														
        ylab("Hg Level (in Liver)") + xlab("Weight (kg)")  

# ----- plot Hg_Liver vs Length
Hg_Liv_weight <- ggplot(peng.dat, aes(x=Length, y=Hg_L)) + 
        stat_smooth(method = "glm", family=Gamma(link="log"), fill = "green", size = 2)+
        geom_point(shape=21, fill="red", size=2, color="#215A89") +
        theme(panel.border = element_rect(fill=NA,color="black", size=2.0,linetype="solid"),legend.position="none") + 				
        theme(panel.background = element_rect(fill='white', colour='gray')) +														
        ylab("Hg Level (in Liver)") + xlab("Length (cm)")  


# ----- Plot and save all combined 
#png(filename="/Users/jfu/Documents/Projects & Manuscripts/manuscripts/Mercury Penguins/GLM_plots_ggplot2.png", width=8.5, height=7, res=1500, units="in")
grid.arrange(Hg_Mus_weight, Hg_Mus_length,Hg_Liv_length, Hg_Liv_weight, ncol = 2)
#dev.off()

# ----- Scatter Plot Hg in Liver and Muscle for each YEAR -------- ----- 
p2006 <- subset(peng.dat, Year == "2006"); p2008 <- subset(peng.dat, Year == "2008"); p2012 <- subset(peng.dat, Year == "2012")
HgM_2006 <- p2006$Hg_M;   HgM_2008 <- p2008$Hg_M;   HgM_2012 <- p2012$Hg_M
 shapiro.test(HgM_2006) ;     shapiro.test( HgM_2008);     shapiro.test(HgM_2012)
 
tiff(filename="/Users/jfu/Documents/Projects & Manuscripts/manuscripts/Mercury Penguins/Plots_TL & Weight_Year_peng2.tiff", width=8, height=7, res=1500, units="in")
dev.off()
      par(mfrow=c(2,6))
      
      par(mfrow=c(2,6), mar=c(3.5,2.5,1,0.5),mgp=c(1.6,0.7,0), oma=c(10,2,3,1))
      
   #Hg in Muscle    --   No clear reletionship
        plot(p2006$Hg_M~p2006$Weight, xlab="Weight (Kg)", ylab="Mercury in Muscle", main="Hg_M vs Weight (2006)");  abline(lm(p2006$Hg_M~ p2006$Weight), col="red", lwd=2);  lines(lowess(p2006$Weight, p2006$Hg_M), col="blue") 
        plot(p2008$Hg_M~p2008$Weight, xlab="Weight (Kg)", ylab="Mercury in Muscle", main="Hg_M vs Weight (2008)"); abline(lm(p2008$Hg_M~ p2008$Weight), col="red", lwd=2);  lines(lowess(p2008$Weight, p2008$Hg_M), col="blue")
        plot(p2012$Hg_M~p2012$Weight, xlab="Weight (Kg)", ylab="Mercury in Muscle", main="Hg_M vs Weight (2012)"); abline(lm(p2012$Hg_M~ p2012$Weight), col="red", lwd=2);  lines(lowess(p2012$Weight, p2012$Hg_M), col="blue")
        plot(p2006$Hg_M~p2006$Length, xlab="Length (cm)", ylab="Mercury in Muscle", main="Hg_M vs Length (2006)"); abline(lm(p2006$Hg_M~ p2006$Length), col="red", lwd=2);  lines(lowess(p2006$Length, p2006$Hg_M), col="blue")
        plot(p2008$Hg_M~p2008$Length, xlab="Length (cm)", ylab="Mercury in Muscle", main="Hg_M vs Length (2008)"); abline(lm(p2008$Hg_M~ p2008$Length), col="red", lwd=2);  lines(lowess(p2008$Length, p2008$Hg_M), col="blue")
        plot(p2012$Hg_M~p2012$Length, xlab="Length (cm)", ylab="Mercury in Muscle", main="Hg_M vs Length (2012)"); abline(lm(p2012$Hg_M~ p2012$Length), col="red", lwd=2);  lines(lowess(p2012$Length, p2012$Hg_M), col="blue")
  
   #Hg in Liver    -   No clear reletionship
       plot(p2006$Hg_L~p2006$Weight, xlab="Weight (Kg)", ylab="Mercury in Liver", main="Hg_L vs Weight (2006)");  abline(lm(p2006$Hg_L~ p2006$Weight), col="red", lwd=2);  lines(lowess(p2006$Weight, p2006$Hg_L), col="blue") 
       plot(p2008$Hg_L~p2008$Weight, xlab="Weight (Kg)", ylab="Mercury in Liver", main="Hg_L vs Weight (2008)");  abline(lm(p2008$Hg_L~ p2008$Weight), col="red", lwd=2);  lines(lowess(p2008$Weight, p2008$Hg_L), col="blue") 
       plot(p2012$Hg_L~p2012$Weight, xlab="Weight (Kg)", ylab="Mercury in Liver", main="Hg_L vs Weight (2012)");  abline(lm(p2012$Hg_L~ p2012$Weight), col="red", lwd=2);  lines(lowess(p2012$Weight, p2012$Hg_L), col="blue") 
       plot(p2006$Hg_L~p2006$Length, xlab="Length (cm)", ylab="Mercury in Liver", main="Hg_L vs Length (2006)");  abline(lm(p2006$Hg_L~ p2006$Length), col="red", lwd=2);  lines(lowess(p2006$Length, p2006$Hg_L), col="blue")     
       plot(p2008$Hg_L~p2008$Length, xlab="Length (cm)", ylab="Mercury in Liver", main="Hg_L vs Length (2008)");  abline(lm(p2008$Hg_L~ p2008$Length), col="red", lwd=2);  lines(lowess(p2008$Length, p2008$Hg_L), col="blue") 
       plot(p2012$Hg_L~p2012$Length, xlab="Length (cm)", ylab="Mercury in Liver", main="Hg_L vs Length (2012)");  abline(lm(p2012$Hg_L~ p2012$Length), col="red", lwd=2);  lines(lowess(p2012$Length, p2012$Hg_L), col="blue") 
       
      
       

# ------ BoxPlot Hg and Biometry across Sex ---------------------- -------
       penguin_sex2 <- filter(peng.dat, Sex != "ND")
       
       Sex.HgM <- ggplot(data = penguin_sex2, aes(x = Sex, y = Hg_M)) +
               geom_boxplot(aes(fill=factor(Sex)),fill="#66A1D2", colour="#27258D", alpha=0.45) +
               theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_blank(), axis.title.x = element_blank()) +
               ylab("Hg Level (in Muscle)")+ xlab("Sex") 	+
               theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  
       # + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)
       
       Sex.HgL <- ggplot(data = penguin_sex2, aes(x = Sex, y = Hg_L)) +
               geom_boxplot(aes(fill=factor(Sex)),fill="#66A1D2", colour="#27258D", alpha=0.45) +
               theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_blank(), axis.title.x = element_blank()) +
               ylab("Hg Level (in Liver)")+ xlab("Sex") 	+
               theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  
       # + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)
       
       Sex.Length <- ggplot(data = penguin_sex2, aes(x = Sex, y = Length)) +
               geom_boxplot(aes(fill=factor(Sex)),fill="#66A1D2", colour="#27258D", alpha=0.45) +
               theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_blank(), axis.title.x = element_blank()) +
               ylab("Length (cm)")+ xlab("Sex") 	+
               theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  
       # + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)
       
       Sex.Weight <- ggplot(data = penguin_sex2, aes(x = Sex, y = Weight)) +
               geom_boxplot(aes(fill=factor(Sex)),fill="#66A1D2", colour="#27258D", alpha=0.45) +
               theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_blank(), axis.title.x = element_blank()) +
               ylab("Weight (kg)")+ xlab("Sex") 	+
               theme(panel.border = element_rect(fill=NA,color="black", size= 1.5, linetype="solid"),legend.position="none")  
       # + geom_point(shape=21, fill="#60D6A7", size=1, color="#215A89", alpha=0.35)
       
svg(filename="/Users/jfu/Documents/Projects & Manuscripts/manuscripts/Mercury Penguins/boxplot_ggplot_SEX.svg", width=8.5, height=7)
grid.arrange(Sex.HgM, Sex.Length, Sex.HgL, Sex.Weight, ncol = 2)
dev.off()

# ---- Map Comparizon Hg Penguins  ------------------------------- ----------

require("maps")
# World map, using geom_path instead of geom_polygon
world_map1 <- map_data("world")
world_sphere1 <- ggplot() +
        theme_bw(20) +  # Here define the width of the the plot box, I think!
        geom_polygon(data = world_map1, aes(x=long, y=lat, group = group), fill = "#E0E0E0", color="#E0E0E0")+
        geom_path(colour = "red") +
        theme(panel.background = element_rect(fill = "#81d4faff", colour = "#81d4faff"),
              #    panel.border = element_rect(size = 2, colour = "black"), 
              panel.grid.major = element_line(colour = "#4fc3f7"),
              panel.grid.minor = element_blank(),
              axis.ticks = element_blank(), 
              axis.title.x = element_text (size = 18, vjust = -0.1), 
              axis.title.y = element_text (size = 18, vjust = 1), 
              axis.text.x = element_text (size = 16, vjust = 0), 
              axis.text.y = element_text (size = 16, hjust = 1.3)) + 
        scale_y_continuous(breaks = (-2:2) * 30) +
        scale_x_continuous(breaks = (-4:4) * 45)

# world_sphere + coord_map("ortho", orientation = c(-49, -9, 0))
# Antarctic coast use this: (-90, -49, 0))
pen_glob <-read.csv("/Users/jfu/Documents/Projects & Manuscripts/manuscripts/Mercury Penguins/Data_&_Rcode/table_compariz.csv",header=T,sep=";") # Loading data
str(pen_glob)

colorpal<-c("#f06292", "#ba68c8", "black", "#64b5f6", "#26c6da", "#81c784", "#dce775", "#ffd54f", "#ff8a65", "#90a4ae")
black <- c("black","black","black","black","black","black","black","black","black","black")
maps1 <- world_sphere1 +
        geom_point(data=pen_glob, 
        aes(x=Long, y=Lat, group=Species, fill=Species, color=Species), 
        size=3.0, shape=21, alpha=.8) +
        scale_color_manual(values=black) + # define color for each symbol according to the categorical values of Species.
        scale_shape_manual(values=c(21,21,21,21,21,21,21,21,21,21)) +  #  symbol shape based on ggplot simbol styles
        scale_fill_manual(values=colorpal)  # # color fill the symbols defined above, if they can be filled (see symbols table for Ggplot2) 


#  svg(filename="world_map_Hg_Peng.svg", width=14, height=6)
# maps1
# dev.off()

Sulax <- ggplot(data=pen_glob, aes(x=Code, y=plot_values, group=Species, fill=Species, color=Species)) +
        geom_bar(position=position_dodge(), stat="identity") +# By default, uses stat="bin", which gives the count in each category
scale_fill_manual(values=colorpal)  +# # color fill the symbols defined above, if they can be filled (see symbols table for Ggplot2) 
scale_color_manual(values=black)  # define color for each symbol according to the categorical values of Species.
        

#  svg(filename="bar_plot_comp.svg", width=14, height=12)
# Sulax 
# dev.off()
Sulax <- ggplot(data=pen_glob, aes(x=Code, y=plot_values)) +

# ---- Happy End ---\ö/--\ö/--\ö/--\ö/--\ö/--\ö/--\ö/--\ö/--\ö/--- ----------



#ancova example -------------

part <- as.numeric(c(3,2,5,2,2,2,7,2,4,7,5,3,4,4,7,5,4,9,2,6,3,4,4,4,6,4,6,2,8,5))
mate <- as.numeric(c(4,1,5,1,2,2,7,4,5,5,3,1,2,2,6,4,2,1,3,5,4,3,3,2,0,1,3,0,1,0))
dose <- factor(c("PLAC", "PLAC", "PLAC", "PLAC", "PLAC", "PLAC", "PLAC", 
                 "PLAC", "PLAC", "LOW", "LOW", "LOW", "LOW", "LOW", "LOW", 
                 "LOW", "LOW", "HIGH", "HIGH",  "HIGH",  "HIGH",  "HIGH",  
                 "HIGH",  "HIGH",  "HIGH",  "HIGH",  "HIGH",  "HIGH",  "HIGH",  "HIGH"))  
libido <- data.frame(part, mate, dose)
str(libido)

hist(libido$part)
hist(libido$mate)
tapply(part, dose, FUN=mean)
tapply(mate, dose, FUN=mean)
tapply(part, dose, sd)
tapply(mate, dose, sd)
#(largest sd)2 / (smalestsd)2 < 2 to check variance
(2.115)^2 / (1.457)^2 # higher than 2!!! Unequal variance Use also levenes
 
#Assumption for the ANCOVA 
# First. Independence of the covariate over treatment effects 
# Second is homogeneity of regression slopes

summary(aov(libido$mate ~ libido$dose))# first assumption
#to analyse the second assumption I would need to create three scatterplots with the linear regression ine plotted
# I'd need to do one for x=mate, y= part, and placebo
# We can also run an ANOVA of participants libido over dose
summary(aov(libido$part ~ libido$dose))# first assumption
model1 <- aov(part~mate+dose, data = libido)
summary(model1)

posthoc_lib <- glht(model1, linfct=mcp(dose ="Tukey"))
summary(posthoc_lib) #Tukey
confint(posthoc_lib)
plot(model1)

 # We need a type III sum of squares calculation
drop1(model1, ~., test="F")
#There is a significant effect of Viagra on the 
#..level of libido when controlling for the covariate 
#..(level of libido) of the mate ), F(2,26)=4.14, p < .05.
# Verifying adjusted mean controling for the covariates
library("effects")
adj.means.lib <- effect("dose", model1, se=TRUE)
summary(adj.means.lib);  
adj.means.lib$se 
summary.lm(model1)
# Post hoc analysis
in
install.packages('multcomp') #  library('multcomp')
install.packages('HH')  #   library('HH')
install.packages('sm')  #   library('sm')

model1$tall <- with(data, interaction(dose, sep='x'))
m1 <- ancova(part~ dose + mate, data=libido)
summary(glht(m1, linfct=mcp(dose="Tukey")))
#http://stats.stackexchange.com/questions/41270/nonparametric-equivalent-of-ancova-for-continuous-dependent-variables




# ---- # ------ ANCOVA Hg Liver Years Rio -------------

peng_year_rj$body_cond <- peng_year_rj$Weight/peng_year_rj$Length
plot(peng_year_rj$Weight, peng_year_rj$body_cond)

hist(peng_year_rj$Hg_L); hist(peng_year_rj$Length)
tapply(peng_year_rj$Hg_L, peng_year_rj$Year, FUN=mean)
tapply(peng_year_rj$Hg_L, peng_year_rj$Year, sd)
tapply(peng_year_rj$Length, peng_year_rj$Year, FUN=mean)
tapply(peng_year_rj$Length, peng_year_rj$Year, sd)
tapply(peng_year_rj$Weight, peng_year_rj$Year, FUN=mean)
tapply(peng_year_rj$Weight, peng_year_rj$Year, sd)

#(largest sd)2 / (smalestsd)2 < 2 to check variance
(2.359)^2 / (0.322)^2 # higher than 2!!! Unequal variance Use also levenes

# LIVER - LIVER - LIVER
head(peng_year_rj)
#peng_year_rj$Hg_Log <- log(peng_year_rj$Hg_L)
model1.hgl.y <- aov(Hg_L ~ Year + Length + Weight , data= peng_year_rj)
drop1(model1.hgl.y, ~., test="F")
plot(peng_year_rj$Hg_L ~peng_year_rj$Year)
text(1, 1.53, "o", cex = .8, col="blue") #x,y
text(2, 2.54, "o", cex = .8, col="blue") #x,y
text(3, 0.70, "o", cex = .8, col="blue") #x,y
text(1, 1.99, "+", cex = .8, col="red") #x,y
text(2, 3.09, "+", cex = .8, col="red") #x,y
text(3, 1.24, "+", cex = .8, col="red") #x,y
text(1, 1.06, "+", cex = .8, col="red") #x,y
text(2, 1.98, "+", cex = .8, col="red") #x,y
text(3, 0.15, "+", cex = .8, col="red") #x,y
plot(lm(Length ~Hg_L, peng_year_rj))
plot(Hg_M ~Length, peng_year_rj)

plot(peng_year_rj$Weight ~peng_year_rj$Year, outline=FALSE)
text(1, 2, "+", cex = .8) #x,y
# Verifying adjusted mean controling for the covariates
library("effects")
adj.means.hgl.y <- effect("Year", model1.hgl.y, se=TRUE)
summary(adj.means.hgl.y);  adj.means.hgl.y$se 

posthoc.hgl.y <- glht(model1.hgl.y, linfct=mcp(Year ="Tukey"))
summary(posthoc.hgl.y) #Tukey
confint(posthoc.hgl.y)
plot(posthoc.hgl.y)

# Posthoc test
summary(glht(model1.hgl.y, linfct=mcp(Year="Tukey")))
model.anov <- aov(model1.hgl.y$Hg_Llog ~ peng_year_rj$Year)

# MUSCLE - MUSCLE - MUSCLE
hist(peng_year_rj$Hg_M); hist(peng_year_rj$Length)
tapply(peng_year_rj$Hg_M, peng_year_rj$Year, FUN=mean)
tapply(peng_year_rj$Hg_M, peng_year_rj$Year, sd)
plot(peng_year_rj$Hg_M ~peng_year_rj$Year)
text(1, 0.60, "o", cex = .8, col="blue") #x,y
text(2, 0.80, "o", cex = .8, col="blue") #x,y
text(3, 0.57, "o", cex = .8, col="blue") #x,y
text(1, 0.65, "+", cex = .8, col="red") #x,y
text(2, 0.87, "+", cex = .8, col="red") #x,y
text(3, 0.63, "+", cex = .8, col="red") #x,y
text(1, 0.54, "+", cex = .8, col="red") #x,y
text(2, 0.74, "+", cex = .8, col="red") #x,y
text(3, 0.51, "+", cex = .8, col="red") #x,y


model1.hgm.y <- aov(sqrt(Hg_M)~ Year  + Length +Weight , data= peng_year_rj)
drop1(model1.hgm.y, ~., test="F")
library("effects")
adj.means.hgm.y <- effect("Year", model1.hgm.y, se=TRUE)
summary(adj.means.hgm.y);  adj.means.hgm.y$se 
posthoc.hgm.y <- glht(model1.hgm.y, linfct=mcp(Year ="Tukey"))
summary(posthoc.hgm.y) #Tukey
confint(posthoc.hgm.y)
plot(posthoc.hgm.y)



# ---- # ------ ANCOVA Hg Liver LATITUDE for 2008 -------------
peng_area_2008$body_cond <- peng_area_2008$Weight/peng_area_2008$Length

# LIVER - LIVER - LIVER
hist(peng_area_2008$Hg_L)
tapply(peng_area_2008$Hg_L, peng_area_2008$Area, FUN=mean)
tapply(peng_area_2008$Hg_L, peng_area_2008$Area, sd)
tapply(peng_area_2008$Length, peng_area_2008$Area, FUN=mean)
tapply(peng_area_2008$Length, peng_area_2008$Area, sd)
tapply(peng_area_2008$Weight, peng_area_2008$Area, FUN=mean)
tapply(peng_area_2008$Weight, peng_area_2008$Area, sd)

plot(peng_area_2008$Weight, peng_area_2008$body_cond)
head(peng_area_2008)
#peng_area_2008$Hg_M <- log(peng_area_2008$HgM)
model1.hgl.a <- aov(Hg_L ~ Area + Length + Weight  , data= peng_area_2008)
drop1(model1.hgl.a, ~., test="F")
plot(peng_area_2008$Hg_L ~peng_area_2008$Area)
text(1, 3.15, "o", cex = .8, col="blue") #x,y
text(2, 2.78, "o", cex = .8, col="blue") #x,y
text(3, 1.85, "o", cex = .8, col="blue") #x,y
text(1, 4.12, "+", cex = .8, col="red") #x,y
text(2, 3.53, "+", cex = .8, col="red") #x,y
text(3, 2.94, "+", cex = .8, col="red") #x,y
text(1, 2.18, "+", cex = .8, col="red") #x,y
text(2, 2.04, "+", cex = .8, col="red") #x,y
text(3, 0.76, "+", cex = .8, col="red") #x,y
# Verifying adjusted mean controling for the covariates
library("effects")
adj.means.hgl.a <- effect("Area", model1.hgl.a, se=TRUE)
summary(adj.means.hgl.a);  adj.means.hgl.a$se 

# Posthoc test
posthoc.hgl.a <- glht(model1.hgl.a, linfct=mcp(Area ="Tukey"))
summary(posthoc.hgl.a) #Tukey
confint(posthoc.hgl.a)
plot(posthoc.hgl.a)

# MUSCLE - MUSCLE - MUSCLE
head(peng_area_2008)
model1.hgm.a <- aov(Hg_M ~ Area  + Length + Weight , data= peng_area_2008)
drop1(model1.hgm.a, ~., test="F")
library("effects")
adj.means.hgm.a <- effect("Area", model1.hgm.a, se=TRUE)
summary(adj.means.hgm.a);  adj.means.hgm.a$se 
plot(peng_area_2008$Hg_M ~peng_area_2008$Area)
text(1, 0.44, "o", cex = .8, col="blue") #x,y
text(2, 0.74, "o", cex = .8, col="blue") #x,y
text(3, 0.54, "o", cex = .8, col="blue") #x,y
text(1, 0.60, "+", cex = .8, col="red") #x,y
text(2, 0.86, "+", cex = .8, col="red") #x,y
text(3, 0.72, "+", cex = .8, col="red") #x,y
text(1, 0.28, "+", cex = .8, col="red") #x,y
text(2, 0.62, "+", cex = .8, col="red") #x,y
text(3, 0.36, "+", cex = .8, col="red") #x,y
# Posthoc test
posthoc.hgm.a <- glht(model1.hgm.a, linfct=mcp(Area ="Tukey"))
summary(posthoc.hgm.a) #Tukey
confint(posthoc.hgm.a)
plot(posthoc.hgm.a)





xyplot(Hg_Llog ~ Length, data=peng_year_rj, groups=Year, aspect="iso", type=c("p","r"),
       auto.key=list(space="right", lines=TRUE, points=FALSE))

ancova.np <- sm.ancova(peng_year_rj$Hg_Llog, peng_year_rj$Length, v$Year, model="equal")
TukeyHSD(ancova.np)
summary(ancova.np)
help()

ancova.np$model


#  Statistics HELL --------
viag.dat <- read.delim("/Users/jfu/Desktop/ViagraCovariate.dat", header=TRUE)
str(viag.dat)
dose <- c(rep(1,9), rep(2,8), rep(3,13))
viag.dat$dose <- factor(dose, levels=c(1:3), labels=c("Placebo","Low Dose","High Dose"))

boxplot(libido~dose, data=viag.dat)
boxplot(partnerLibido~dose, data=viag.dat)
viag.mod1 <- aov(libido~dose, data=viag.dat)
summary(viag.mod1)
viag.mod2 <- aov(libido~dose+partnerLibido, data=viag.dat)
summary(viag.mod2)
drop1(viag.mod2,~.,test="F") # type III SS and F Tests
#The peson's libido is influenced by their partner's libido. 
#..When the effect of partner's libido is removed, the effect of 
#.. Viagra is significant (p 0.027).

library("effects")
adjusted.means <- effect("dose", viag.mod2, se=TRUE)
summary(adjusted.means)
adjusted.means$se 

summary.lm(viag.mod2)
posthocs <- glht(viag.mod2, linfct=mcp(dose="Tukey"))
summary(posthocs)
confint(posthocs)

model = lm(peng_area_2008$Hg_L  ~ log(peng_area_2008$Length) )
summary(model)
plot(log(peng_area_2008$Length), peng_area_2008$Hg_L) #peng_year_rj
abline(model)
cor(peng_year_rj$Hg_M, peng_year_rj$Weight)
















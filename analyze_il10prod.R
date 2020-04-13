dat <- read.csv("il10prod.csv")#,stringsAsFactors = FALSE)
basey <- dat[dat$cond=="ref",]$y
dat$y <- dat$y - basey
dat <- dat[-1,]

dat$y[dat$replicate==1] <- dat$y[dat$replicate==1]*0.02/(893-basey)
dat$y[dat$replicate==2] <- dat$y[dat$replicate==2]*0.05/(893-basey)
#893 = 2%, in plus
#893 = 5%, in minus

var.test(
  dat$y[dat$cond=="naive" & dat$isko=="YES"],
  dat$y[dat$cond=="naive" & dat$isko=="NO"],
  alternative="greater"
)

var.test(
  dat$y[dat$cond=="nb" & dat$isko=="YES"],
  dat$y[dat$cond=="nb" & dat$isko=="NO"],
  alternative="greater"
)




datnb <- dat[dat$cond=="nb",]
linmod <- lm(y~isko+foxp3,datnb)
anova(linmod)

dat
datnb <- dat[dat$cond=="nb" & dat$foxp3=="plus",]
t.test(
  datnb$y[datnb$isko=="YES"],
  datnb$y[datnb$isko=="NO"])

datnb <- dat[dat$cond=="nb" & dat$foxp3=="minus",]
t.test(
  datnb$y[datnb$isko=="YES"],
  datnb$y[datnb$isko=="NO"])



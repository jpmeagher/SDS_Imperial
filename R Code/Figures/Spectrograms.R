## SDS: Figure 1
library(ancestralbat)
library(signal)

# Get example spectrograms
w <- 512
anpa <- specgram(mexicancalls$chirp1[[2]], w, Fs = 500000, window = hamming(w),
                 overlap = 15*(w/16))

bapl <- specgram(mexicancalls$chirp1[[72]], w, Fs = 500000, window = hamming(w),
                 overlap = 15*(w/16))

myyu <- specgram(mexicancalls$chirp1[[248]], w, Fs = 500000, window = hamming(w),
                 overlap = 15*(w/16))

ptpa <- specgram(mexicancalls$chirp1[[366]], w, Fs = 500000, window = hamming(w),
                 overlap = 15*(w/16))

# Produce Image for each and save as pdf
((23 - max(bapl$t*1000) )/ 2)

x11()
anpa_s <- 10*log10(abs(t(anpa$S^2)))
anpa_s[anpa_s < 0] <- 0
image(anpa$t*1000, anpa$f/1000, anpa_s, col = grey(512:0 /512), ylim = c(9, 212), 
      xlim = c(0 - ((23 - max(anpa$t*1000) )/ 2), max(anpa$t*1000) + ((23 - max(anpa$t*1000) )/ 2)),
      xlab = "Time (ms)", ylab = "Frequency (kHz)", cex.lab = 1.5, cex.axis = 1.5)
box()

x11()
bapl_s <- 10*log10(abs(t(bapl$S^2)))
bapl_s[bapl_s < 0] <- 0
image(bapl$t*1000, bapl$f/1000, bapl_s, col = grey(512:0 /512), ylim = c(9, 212), 
      xlim = c(0 - ((23 - max(bapl$t*1000) )/ 2), max(bapl$t*1000) + ((23 - max(bapl$t*1000) )/ 2)),
      xlab = "Time (ms)", ylab = "Frequency (kHz)", cex.lab = 1.5, cex.axis = 1.5)
box()

x11()
myyu_s <- 10*log10(abs(t(myyu$S^2)))
myyu_s[myyu_s < 0] <- 0
image(myyu$t*1000, myyu$f/1000, myyu_s, col = grey(512:0 /512), ylim = c(9, 212), 
      xlim = c(0 - ((23 - max(myyu$t*1000) )/ 2), max(myyu$t*1000) + ((23 - max(myyu$t*1000) )/ 2)),
      xlab = "Time (ms)", ylab = "Frequency (kHz)",cex.lab = 1.5, cex.axis = 1.5)
box()

x11()
ptpa_s <- 10*log10(abs(t(ptpa$S^2)))
ptpa_s[ptpa_s < 0] <- 0
image(ptpa$t*1000, ptpa$f/1000, ptpa_s, col = grey(512:0 /512), ylim = c(9, 212), 
      xlim = c(0 - ((23 - max(ptpa$t*1000) )/ 2), max(ptpa$t*1000) + ((23 - max(ptpa$t*1000) )/ 2)),
      xlab = "Time (ms)", ylab = "Frequency (kHz)", cex.lab = 1.5, cex.axis = 1.5)
box()

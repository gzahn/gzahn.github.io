library(qrcode)
url <- "https://gzahn.github.io"

qr <- qr_code(url,ecl = "H")
plot(qr)

png("~/Desktop/GIT_REPOSITORIES/gzahn.github.io/media/website_qrcode.png",width = 2,height = 2,
    units = "in",res=300)
qrcode::add_logo(code = qr,logo = "~/Desktop/GIT_REPOSITORIES/gzahn.github.io/media/mushroom _clipart.png",ecl = "L") %>% plot
dev.off()

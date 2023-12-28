library(hexSticker)
imgurl <- "~/Documents/hub_taxa.png"

sticker(imgurl, package = "hubfindr",filename = "~/Desktop/hex.png",
        s_width = .5,
        s_height = .5,
        s_x = 1,
        s_y = .8,
        p_y = 1.5,
        p_size = 16,
        h_fill = "#173021",
        spotlight = TRUE,
        l_y = .8,
        l_alpha = .25,
        l_width = 100,
        l_height = 5)
getwd()

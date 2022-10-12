library(hexSticker)

imgurl <- "./media/stickergraph.png"

sticker(imgurl, package = "Microbiome Boot Camp",filename = "./media/hex.png",
        s_width = .8,
        s_height = .5,
        s_x = .9,
        s_y = .8,
        p_y = 1.4,
        p_size = 11,
        h_fill = "black",
        spotlight = TRUE,
        l_y = 1,
        l_alpha = .25,
        l_width = 100,
        l_height = 5)

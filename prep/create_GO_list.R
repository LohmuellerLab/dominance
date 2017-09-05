library(GO.db)

mf <- get("GO:0003674", GOMFOFFSPRING)
cc <- get("GO:0005575", GOCCOFFSPRING)
bp <- get("GO:0008150", GOBPOFFSPRING)

binding <- get("GO:0005488", GOMFOFFSPRING)
catalytic <- get("GO:0003824", GOMFOFFSPRING)
structural <- get("GO:0005198", GOMFOFFSPRING)

write(mf, "data/reference/MF.txt")
write(cc, "data/reference/CC.txt")
write(bp, "data/reference/BP.txt")

write(binding, "data/reference/binding.txt")
write(catalytic, "data/reference/catalytic.txt")
write(structural, "data/reference/structural.txt")
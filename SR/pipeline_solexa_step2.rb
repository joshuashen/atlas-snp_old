## 
prefix = ARGV[0]  # prefix of output 

ref = ARGV[1]  # reference .bfa file

maps = ARGV[2..-1].join("\t")   ## all the .map files, DIR/s*.map would do

# puts maps

# merge all .map files
system("maq mapmerge #{prefix}.map #{maps}")

# assemble
system("maq assemble #{prefix}.cns #{ref} #{prefix}.map 2> #{prefix}.cns.log")

# call SNP
system("maq cns2snp #{prefix}.cns > #{prefix}.cns.SNP")

# pileup
system("maq pileup -sP #{ref} #{prefix}.map > #{prefix}.pileup")

# mapview
system("maq mapview -bN #{prefix}.map > #{prefix}.map.mapview")


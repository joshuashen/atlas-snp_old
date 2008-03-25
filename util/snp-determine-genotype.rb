## call genotype after determine het_hom

while line=ARGF.gets do 
  cols=line.split(/\s+/)
  
  if cols[9].to_f >= cols[5].to_f*0.9  # or   cols[-1].to_f <0.01 # homo
    gt = "#{cols[6]}#{cols[6]}"
  else
    # het
    ref = cols[2].upcase
    gt = [ref, cols[6]].sort.join()
  end

  puts "#{line.chomp}\t#{gt}"
end


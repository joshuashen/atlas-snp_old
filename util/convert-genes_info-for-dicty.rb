## 


$chromosome_change={"M"=>"DDB0169550","2F"=>"DDB0215018", "BF"=>"DDB0220052", "3F"=>"DDB0215151", "1"=>"DDB0232428", "2"=>"DDB0232429","3"=>"DDB0232430", "4"=>"DDB0232431", "5"=>"DDB0232432", "6"=>"DDB0232433"}


while line = ARGF.gets do 
  cols = line.split(/\s+/)

  
  chr = cols[2]
  if $chromosome_change.key?(chr)
    puts "#{cols[0]}\t#{$chromosome_change[cols[2]]}\t#{cols[3]}\t#{cols[4]}\t#{cols[1]}\t#{cols[5..-1].join(" ")}"
  end
end

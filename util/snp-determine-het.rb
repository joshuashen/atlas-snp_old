## OLD:  minor allele >20 % --> heterozygous SNP
## NEW: based on probability from binomial distribution

# Null hypothesis: the SNP is heterozygous, alterative: the SNP is homozygous. 
## suggest Cutoff value: if Pr < 0.02, reject the null hypothesis ==> the SNP is homozygous.

def factorial(k)
  if k <= 1 
    return 1
  else
    return k * factorial(k-1)
  end
end

def combinatory(n,i)
  t = 1.0
  1.upto(i) do |j|
    t = t* (n-j+1) / j
  end
  return t
end

while line= ARGF.gets do 
  cols = line.split(/\s+/)
  cov = cols[5].to_f
  alter = cols[9].to_i

  
  ## Binomial distribution:
## F(x; n, 0.5) = Pr(X<= x) = SUM (take j from n)*0.5^n
  pr = 0
  0.upto(alter) do |i|
#    pr += (factorial(cov) / (factorial(i)* factorial(cov-i))) * (0.5**cov )
    pr += combinatory(cov,i) * (0.5**cov)
  end
  
#  if pr < 0.5
#    pr = 1- pr
#  end
  
  pr = 1 - pr
# 1 - pr # probability of error to reject null hypothesis -- that this is a heterozygous SNP
  pr = (pr*10000).round/10000.0
  puts "#{line.chomp}\t#{pr}"

end

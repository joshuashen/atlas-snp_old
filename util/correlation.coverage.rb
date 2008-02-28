## compute correlation coefficient between two lists of genes

# pseudo-code:
# sum_sq_x = 0
# sum_sq_y = 0
# sum_coproduct = 0
# mean_x = x[1]
# mean_y = y[1]
# for i in 2 to N:
#    sweep = (i - 1.0) / i
#    delta_x = x[i] - mean_x
#    delta_y = y[i] - mean_y
#    sum_sq_x += delta_x * delta_x * sweep
#    sum_sq_y += delta_y * delta_y * sweep
#    sum_coproduct += delta_x * delta_y * sweep
#    mean_x += delta_x / i
#    mean_y += delta_y / i 
# pop_sd_x = sqrt( sum_sq_x / N )
# pop_sd_y = sqrt( sum_sq_y / N )
# cov_x_y = sum_coproduct / N
# correlation = cov_x_y / (pop_sd_x * pop_sd_y)

def computeCor(a1, a2)
  s1 = a1.size 
  
  return -100 if s1 < 2 or s1 != a2.size
  
  sumsqx, sumsqy, sumcproduct = 0,0,0
  meanx, meany = a1[0], a2[0]
  1.upto(s1-1) do |i|
    j = i + 1.0
    sweep = (j -1.0) / j
    deltax = a1[i] - meanx
    deltay = a2[i] - meany
    sumsqx += deltax * deltax * sweep
    sumsqy += deltay * deltay * sweep
    sumcproduct += deltax * deltay * sweep
    meanx += deltax / j
    meany += deltay / j
  end
  
  popsdx = (sumsqx / s1.to_f)**0.5
  popsdy = (sumsqy / s1.to_f)**0.5
  covxy = sumcproduct/s1.to_f
  cc = covxy / (popsdx * popsdy)
  return cc
end

list1 = ARGV[0]
list2 = ARGV[1]

genes1 = {}
genes2 = {}

File.new(list1, 'r').each do |line|
  cols = line.split(/\s+/)
  gene, c = cols[0], cols[2].to_f
  if !genes1.key?(gene)
    genes1[gene] = []
  end
  genes1[gene] << c
end

File.new(list2, 'r').each do |line|
  cols = line.split(/\s+/)
  gene,c = cols[0], cols[2].to_i
  if !genes2.key?(gene)
    genes2[gene] = []
  end
  genes2[gene] << c
end

genes1.keys.sort.each do |gene|
  cor = computeCor(genes1[gene], genes2[gene])
  puts "#{gene}\t#{cor}\t#{genes1[gene].size}"
end


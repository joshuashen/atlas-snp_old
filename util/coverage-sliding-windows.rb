## 

# input format: 

# ENm010  1955    9       chr7    27126000

def mean(a)
  return -1 if a.size == 0
  t = 0
  a.map {|i| t+= i}
  
  return t/a.size.to_f
end

pqueue = []
cqueue = []

while line = ARGF.gets do 
  cols = line.split(/\s+/)
  
  name, pos, cov, chr, gpos = cols[0], cols[1].to_i, cols[2].to_i, cols[3], cols[4].to_i
  
  pqueue << gpos
  cqueue << cov

  if pqueue.size >= 1000
    pmean = mean(pqueue)
    cmean = mean(cqueue)
    if pmean > 0
      puts "#{name}\t#{chr}\t#{pmean}\t#{cmean}"
    end
    pqueue = []
    cqueue = []
  end
  
end
  
pmean = mean(pqueue)
cmean = mean(cqueue)
if pmean >=0 
  puts "#{name}\t#{chr}\t#{pmean}\t#{cmean}"
end

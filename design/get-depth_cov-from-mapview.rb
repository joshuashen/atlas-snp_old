# purpose: estimate the relationship between shape and avg depth-coverage. 
## shape is the second parameter from negative binomial distribution (over-dispersed Poisson)

$readsize  = 35  ## determine the window size
$minQual = 10 # min mapping quality
$chrSize = {}

def main 
  mapview=ARGV[0]
  refstat = ARGV[1]  # size of chromosomes

  getSize(refstat)
  readMapview(mapview)
end

def getSize(refstat)

  File.new(refstat, 'r').each do |line|
    if line=~ /^(\S+)\s+(\d+)/
      chr, s = $1, $2.to_i
      $chrSize[chr] = s

    end
  end
end

def printcov(depcov)
  return if depcov.size < 1
  puts depcov.join("\n")
end

def readMapview(mapview)
  current = 0
  prevchr = ''
  depCov = []

  File.new(mapview, 'r').each do |line|
    cols = line.split(/\s+/)
    chr,pos,flag,qual = cols[1], cols[2].to_i, cols[5], cols[6].to_i
    
    # only take good pairs
    next if flag != "18" or qual < $minQual or !($chrSize.key?(chr))
    
    if chr != prevchr # a new chr
      $stderr.puts "#{chr}"
      printcov(depCov)
      depCov = Array.new($chrSize[chr],0)
      
    end
    prevchr = chr
    
    pos.upto(pos + 35 - 1) do |i|
      depCov[i-1] += 1
    end
    
  end

  printcov(depCov)

end

main()

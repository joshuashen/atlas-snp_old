## call genotypes at SNP sites based on Maximum Likelihood. 
# the paramemter in question is the genotype.

require 'getoptlong'

opts = GetoptLong.new(
    ["--snp", "-s", GetoptLong::REQUIRED_ARGUMENT],
    ["--hapmap", "-m", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help")
  $stderr.puts "Usage: ruby __.rb -s snp.list [ -m hapmap.maf ] "
  $stderr.puts " -s snp.list, the output of atlas-snp.rb"
  $stderr.puts " -m hapmap.maf, HapMap-reported SNP allele frequency"
  exit
end

def likelihood(refbase,snpbase,refbasenum,snpbasenum, otherbasenum, chr, position, logisticOdds)
  # only consider three states: AA, AB, BB, A as the reference allele
  # compute llh on condition of binomial p = 0.4, 0.5, and 0.6 

  # example:  AAAT, A as the ref, T's frequency is 0.2 in the population, 
  # and logistic  l_p = 0.2, assuming p=0.5, and overall error rate is 0.1%,
  #  then:
  # llh(AAAT|AT) = (1-0.001)^4*0.5^4 = 0.0625 based on binomial 
  # llh(AAAT|AA) = (1-0.001)^3 * 0.001/(0.001 + 2*0.2) = 0.0025   
  #     (based on logistic/Bayesian)
  # llh(AAAT|TT) =  ignore
  # --> Odds_ratio = 25
  # and log10(OR) = 1.4 
  # most likely genotype: AT

  llh = {'AA' => 1, 'AB' => 1, 'BB' => 1}

  # check if the snpbase is of known allele from HapMap:
  if $population.key?(chr) and $population[chr].key?(position) and $population[chr][position].key?(snpbase)  
    prior = $population[chr][position][snpbase]
  else
    prior = $priorSNPRate
  end

  # related to ref -> alter prob
  #   logisticOdds = from argument
  # ref -> SNP error prob
  emitAlleleErrorProb = $priorError/($priorError+logisticOdds * prior)
  
  # non-ref_to_SNP error prob
  hiddenlAlleleErrorProb = $priorError/($priorError+logisticOdds * $priorSNPRate)
  
  llh = {}
  # SNPbase are caused by sequencing errors:
  llh['aa'] = (1 - hiddenAlleleErrorProb)**refbasenum * (emitAlleleErrorProb)**snpbasenum
#   llh['ab_4'] = ((1 - hiddenAlleleErrorProb) * 0.4)**refbasenum * (0.6*(1 - hiddenAlleleErrorProb))**snpbasenum
  llh['ab'] = ((1 - hiddenAlleleErrorProb) * 0.5)**(refbasenum+snpbasenum)
#  llh['ab_6'] = ((1 - hiddenAlleleErrorProb) * 0.6)**refbasenum * (0.4*(1 - hiddenAlleleErrorProb))**snpbasenum

  llh['bb'] = (hiddenAlleleErrorProb)**refbasenum * (1-hiddenAlleleErrorProb)**snpbasenum
  
 
  array = llh.keys.sort {|a,b| llh[b] <=> llh[a]}
  lod = Math.log10(llh[array[0]]/(llh[array[1]]+0.0000001*llh[array[0]]))
  if array[0] == 'aa'
    return "#{refbase}#{refbase}\t#{lod}\tref-homo"
  elsif array[0] == 'ab'
    return "#{refbase}#{snpbase}\t#{lod}\thet"
  elsif array[0] == 'bb'
    return "#{snpbase}#{snpbase}\t#{lod}\thomo"
  end

end
 
def logistic
 # compute logistic odds: p(SNP|l_p) / p(error|l_p)
end


File.new(optHash["--snp"], 'r').each do |line|
  cols=line.split(/\s+/)
  next if cols.size < 10
  chr, pos, refbase,homo,cov,snpbase,num,other,readinfo =cols[0], cols[1].to_i, cols[2].upcase, cols[3].to_i,cols[5].to_i, cols[6], cols[9].to_i, cols[10].to_i - cols[9].to_i, cols[11]
  
  print cols[0..-2].join("\t") + "\t"

  if snpbase =="SNPBase"
    puts readinfo + "\tless_strand\tPr_e_c\tstrandP\tGtype\tLOD"
  else
    eva = 0
    eP = 1
    dirs = {'+' => 0, '-' => 0}
    readinfo.split(";").each do |r|
      if r=~ /^(\S+)\((\d+)\)\S+\((\d+)\)\S+([+|-])\S+\((\S+)\/(\S+)\/(\d+)\)(\S+)/
        base,qual,dist,dir,sub,indel,tail,swap = $1,$2.to_i,$3.to_i,$4,$5.to_f,$6.to_f,$7.to_i,$8
        next if base != snpbase
        dirs[dir] += 1
        if base == 'N'
          logOdd = 0
          errorPosteriorFull = 1
        else
          if dist > 120 #  distance to 3' 
            dist = 120 
          end
          
          if swap =="swap"
            swap = 1
          elsif swap =="mnp"
            swap = 0.5
          else
            swap = 0
          end
          
          # apply logistic regression 
          predict = $logistic[:intercept] + $logistic[:swap]*swap + $logistic[:homo]*homo + \
          $logistic[:qual] * qual + $logistic[:dist] * dist
          logisticP = 1 - Math.exp(predict) / (1 + Math.exp(predict))
          # map logisticP to a bin, 0.05, 0.15, ..., 0.95
          bin = (logisticP / 0.1 - 0.5).round 
          
          # boundary adjust
          if bin > 9 
            bin = 9
          elsif bin < 0
            bin = 0
          end
          
          ePrior = $errorPrior[bin]
          sPrior = $snpPrior[bin]
          
          odds = sPrior/ePrior

          gtype= likelihood(refbase,snpbase,cov - num - other, num, other, chr, pos, odds)
          print ""

          ## note: need to use base quality score of each read.? -- treat refbase as quality 30 
        end
      end
    end
  end
end

          

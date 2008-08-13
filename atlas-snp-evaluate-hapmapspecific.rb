## todo: integrate HapMap allele frequency as the Prior SNP rate for known SNPs


## estimate the probability of error of raw SNP calls on the condition of:
##  1. quality scores
#   2. longest homopolymer size in a 13-bp reference window centered on the substitution base
#   3. distance to 3' end on the read
#   4. "swapped" bases  -- specific to 454 sequencing
#  ideas wanted:   *agreement among mapped reads*

## ------------------------------------------------------------------
## I. logistic linear regression, trained on fly 454 FLX data:
#
#      log(p/(1-p)) = a + b*swap + c*homo + d*qual + e * dist_3p
#
##    trained beta:  
#       (Intercept)      y$swap      y$homo      y$qual     y$dist3 
##       1.37124483 -2.08600342 -0.63242288  0.03681018  0.01027486
##
##  p-value: all < e-16
#
# -------------------------------------------------------------------
## II. Bayesian inference:
#
## p(SNP | logisticP) = p(logisticP | SNP) * P(SNP) / P(logisticP)
# or
#  P(error | logisticP) = p(logisticP | error ) * p(error) / P(logisticP)
#
### where p(error) = error rate, about 0.08%
##        p(logisticP | error) = distribution from logistic model
##        p(logisticP)  = SUM(P(logisticP |error) * P(error) + P(logisticP | SNP) * P(SNP) ) and might be a constant

##  --------------------------------------------------------------
## Priors:

## error:
#         *: substitution errors
#         #: SNPs
#
#  (proportion)
#         ^
#         |
#    0.35-|      #
#         |
#         |
#    0.30-|
#         |                                              *
#         |
#    0.25-|
#         |
#         |           #
#    0.20-| #
#         |
#         |
#    0.15-|                                        *                                     
#         |                                                                                       
#         |                                                                                   
#    0.10-|      *    *    #                                                                    
#         |                               *    *                                               
#         | *              *    #                            
#    0.05-|                     *    *
#         |                          #    
#         |                               #    #    #    #
#       0-|________________________________________________ (logistic-P)
#           |    |    |    |    |    |    |    |    |    |
#          0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95
#
#
def combinatory(n,i)
  t = 1.0
  1.upto(i) do |j|
    t = t* (n-j+1) / j
  end
  return t
end


require 'getoptlong'

opts = GetoptLong.new(
    ["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
    ["--errorRate", "-e", GetoptLong::OPTIONAL_ARGUMENT],
    ["--snpRate", "-s", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT],
    ["--frequency", "-f", GetoptLong::NO_ARGUMENT],
    ["--compare","-c", GetoptLong::NO_ARGUMENT],
    ["--dir","-d", GetoptLong::OPTIONAL_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") or !optHash.key?("--input")
  $stderr.puts "Usage: ruby __.rb -i SNP.list [ -e estimated_substitution_error_rate ] [ -s estimated_SNP_rate ] > output"
  $stderr.puts "       -e   default 0.0008; this is optimized for 454 FLX"
  $stderr.puts "       -s   default 0.001, optimized for human "
  $stderr.puts "       -f   enable a component specifical for human hapmap"
  $stderr.puts "       -c   this option should be used together with option '-f'. If '-f' is appended, once '-c' is appended, prior of P(SNP) and Pr(e|c) produced by using the version without '-f' will be appended to the last two column for comparison purpose."
  $stderr.puts "       -d   this option is required when '-f' is turned on. It must be followed by the directory where human hapmap files(foward strand) are stored"
  exit
end

if optHash.key?("--errorRate")
  $errorRate = optHash["--errorRate"].to_f
else
  $errorRate = 0.0008
end

if optHash.key?("--snpRate")
  $snpRate_default = optHash["--snpRate"].to_f
else
  $snpRate_default = 0.001
end

#if optHash.key?("--dir")
 # $dir=optHash["--dir"].to_s
#end

#$stderr.puts $dir

if optHash.key?("--frequency")
  $byhapmap = true
  if optHash.key?("--dir")
   $hapmapPath=optHash["--dir"].to_s
    $stderr.puts $hapmapPath
  else
    $stderr.puts "Couldn't move forward until a path for hapmap is denoted!"
    exit;
  end
  
  if optHash.key?("--compare")
    $comparePrec=true
  else
    $comparePrec=false
  end
else
  $byhapmap = false
  if optHash.key?("--compare")
    $stderr.puts "No comparison needed!"
    exit;
  end
end


$snpRate=0
## based on trainning data from fly
## TODO: let user choose a set of trained parameters from other data sets.
 
$logistic = {:intercept => 1.37, :swap => -2.09, :homo => -0.632, :qual => 0.0368, :dist =>0.0103}

# $errorPrior= {0.05 => 0.062, 0.15 => 0.099, 0.25 => 0.091, 0.35 => 0.059, 0.45 => 0.053, 0.55=> 0.054, 0.65=> 0.073, 0.75 => 0.08, 0.85 => 0.15, 0.95 => 0.28 }

$errorPrior = [0.062, 0.099, 0.091, 0.059, 0.053, 0.054, 0.073, 0.08, 0.15, 0.28]

# $snpPrior = {0.05=> 0.20, 0.15=> 0.35, 0.25 =>0.21, 0.35 => 0.097, 0.45 => 0.063, 0.55 => 0.035, 0.65=> 0.019, 0.75 => 0.011, 0.85 => 0.0064, 0.095 => 0.0032 }

$snpPrior = [0.20, 0.35, 0.21, 0.097, 0.063, 0.035, 0.019, 0.011, 0.0064, 0.0032]


def harmonicMean(array)
  return nil if array.size < 1
  t = 0.0
  array.each do i
    t += 1/i.to_f
  end

  return 1/t
end

#$hapmapDir='/data/pmrs/zwan/atlas-snp/hapmap/'
#This program has strict requirement on hapmap. In this freeze, we only considered the population CEU.  In next freeze, it will be more user friendly to help users to locate the specifical hapmap they want to use.

$hapmapDir=$hapmapPath #'/data/pmrs/zwan/atlas-snp/hapmap/' is a directory storing human hapmap for population CEU.
$prefix='allele_freqs_'
$post='_CEU_r23a_nr.b36_fwd.txt'
$hapmap= Hash.new {|h,k| h[k]=Hash.new}
$hapmapFile=''

def createHMHash()
    File.new($hapmapFile,'r').each do |line|
        cols=line.split(/\s+/)
        if line=~ /^rs#\s+(\S+)/
        else
            chr, coord, freq=cols[1],cols[2].to_i,cols[-3].to_f
            strand=cols[3]
            if strand=='+'
              $hapmap[chr][coord] = freq
            end
        end
    end
end

chr=''
count=0

File.new(optHash["--input"], "r").each do  |line|
  cols=line.split(/\s+/)

  refbase,homo,cov,snpbase,num,alternum,readinfo =cols[2].upcase, cols[3].to_i,cols[5].to_i, cols[6], cols[9].to_i, cols[10].to_i, cols[11]
  
  print cols[0..-2].join("\t") + "\t"
  
  ################Get frequencies of known SNPs from Hapmap#######################
  hmstatus='unknown'
  
  chrom,pos=cols[0],cols[1].to_i
  if chrom!=chr and count>0
    $stderr.puts chr+" has been done!"
  end
  count+=1
  chr=chrom
  if ($byhapmap)
    $hapmapFile=$hapmapDir+$prefix+chrom+$post
    #$stderr.puts $hapmapFile
  end
  
  if ($byhapmap) and File.exist?($hapmapFile)
    if !$hapmap.key?(chrom)
      createHMHash()
      $stderr.puts "Creating hash for hapmap "+chrom+" done!"
    end
    #$stderr.puts "Creating hash for hapmap "+chrom+" done!"
  end
  
  if ($byhapmap)
    if $hapmap[chrom][pos]!=nil
      hmstatus='known'
    end
    if $hapmap[chrom][pos]!=nil and $hapmap[chrom][pos]!=0
      $snpRate=$hapmap[chrom][pos].to_f
    else
      $snpRate=$snpRate_default
    end
  end
  #####################################################
  if snpbase =="SNPBase"
    if ($byhapmap)
      puts readinfo + "\tless_strand\tPr_e_c\tstrandP\tsnpRate\thapmap status\tdefault SNP rate\told Pr_e_c\n"
    else
      puts readinfo + "\tless_strand\tPr_e_c\tstrandP\n"
    end
  else
    eva = 0
    eva_default=0
    eP = 1
    eP_default=1
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
          
          # posterior, omit the demonimator, because a) it's constant, b) we are interested in log Odds.
          errorPosterior = ePrior * $errorRate
          snpPosterior = sPrior * $snpRate.to_f
          
          
          snpPosterior_default = sPrior * $snpRate_default.to_f
          logOdd_default = ( Math.log(snpPosterior_default / errorPosterior) * 100 ).round / 100.to_f
          errorPosteriorFull_default = 1.0 / (1 + snpPosterior_default/errorPosterior)
          if ($byhapmap)
            logOdd = ( Math.log(snpPosterior / errorPosterior) * 100 ).round/ 100.to_f
            errorPosteriorFull = 1.0 / (1 + snpPosterior/errorPosterior)
          end
          
          #logOdd = ( Math.log(snpPosterior / errorPosterior) * 100 ).round/ 100.to_f
        end
        
        eva_default += logOdd_default
        eP_default *=errorPosteriorFull_default
        
        if ($byhapmap)
          eva += logOdd
          eP *= errorPosteriorFull
        end
        if ($byhapmap)
          print "#{r}(#{logOdd})(#{errorPosteriorFull});"
        else
          print "#{r}(#{logOdd_default})(#{errorPosteriorFull_default});"
        end
        
      end
    end

#    eP = (eP*).round/100000.0
    x = [dirs['+'], dirs['-']].min

## Binomial distribution:
## F(x; n, 0.5) = Pr(X<= x) = SUM (take j from n)*0.5^n
    strandPr = 0
    0.upto(x) do |i|
      strandPr += combinatory(num,i) * (0.5**num)
    end
    
    if $byhapmap==true and $comparePrec==true
      print "\t#{x}\t#{eP}\t#{strandPr}\t#{$snpRate}\t#{hmstatus}\t#{$snpRate_default}\t#{eP_default}\n"
    elsif $byhapmap==true and $comparePrec==false
      print "\t#{x}\t#{eP}\t#{strandPr}\t#{$snpRate}\t#{hmstatus}\t--\t--\n"
    elsif $byhapmap==false
      print "\t#{x}\t#{eP_default}\t#{strandPr}\n"
    end
    
  end
end
$stderr.puts chr+" has been done!"
#$selsnpoutput.close

        






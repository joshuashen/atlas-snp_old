require 'getoptlong'
require 'bigdecimal'

#Program overview
#    'atlasSNPEvaluate.rb' is used to evalute substitutions by looking at the posterior probability of
#  a substitution site to be a true SNP (Posterior(SNP|Sj,Cov)j), which is appended in the last column of the output.
#  
#Baysian's framework
#    Posterior(SNP|Sj,Cov)j = Pr(Sj|SNP,Cov)*Prior(SNP|Cov) / ( Pr(Sj|Error,Cov)*Prior(Error|Cov) + Pr(Sj|SNP, Cov)*Prior(SNP|Cov) )
#    Posterior(SNP|Sj,Cov)j: Posterior probability of a substitution site j to be a true SNP, when signal is Sj and variant at specific variant read coverage,cov.
#    Pr(Sj|SNP,cov): Density distribution of Sj for SNPs
#    Pr(Sj|Error): Density distribution of Sj for Errors
#    Prior(SNP|Cov): Estimated SNP rate based variant read coverage
#    Prior(Error|Cov): Estimated error rate based variant read coverage.
#    Sj: Prior probability of a site j to be a true SNP.
#  
#Command line arguments
#  Usage: ruby atlasSNPEvaluate_V3.0.rb -i [Substitution list] -x [option for XLR data] -e [priorerror] -o [output]
#  -i is a required argument, which must be followed by the substitution list.
#  -x is an optional argument. This select must be appended if the input is XLR data.
#  -e is an optional argument with an associate value, Prior(Error|Cov) when variant coverage greater than 2. The default value is 0.1 if this option is ignored.
#  -o is a required argument with an associate output file name with path.
#
#Example:
#  1. Evaluate on the XLR substitution list 'chr1.xm.SNP.list', with Prior(Error|Cov), 0.4 when variant coverage > 2. The
#     output will be in file chr1/chr1.xm.SNP.list.eva.prior_0.4. 
#         ruby atlasSNPEvaluate_V3.0.rb -i chr1.xm.SNP.list -x -e 0.4 -o chr1/chr1.xm.SNP.list.eva.prior_0.4
#  2. Evaluate on the FLX substitution list 'chr1.xm.SNP.list', with default Prior(Error|Cov), 0.1 when variant coverage > 2.
#     The output will be in file chr1/chr1.xm.SNP.list.eva.prior_0.1.
#         ruby atlasSNPEvaluate_op_V3.1.rb -i chr1.xm.SNP.list -o chr1/chr1.xm.SNP.list.eva.prior_0.1
#
#
#
class AtlasEvaluate
#Initialize options
def AtlasEvaluate.processArguments()
opts = GetoptLong.new(
    ["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT],
    ["--XLR", "-x", GetoptLong::NO_ARGUMENT],
    ["--output","-o", GetoptLong::REQUIRED_ARGUMENT],
    ["--priorerror","-e", GetoptLong::OPTIONAL_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

AtlasEvaluate.usage() if (optHash.key?("--help"));

AtlasEvaluate.usage() if (optHash.empty?);
return optHash
end

#Usage information
def AtlasEvaluate.usage(msg='')
  unless (msg.empty?)
        puts "\n#{msg}\n"
  end
  puts "
Program overview
    'atlasSNPEvaluate.rb' is used to evalute substitutions by looking at the posterior probability of
  a substitution site to be a true SNP (Posterior(SNP|Sj,Cov)j), which is appended in the last column of the output.
  
Baysian's framework
    Posterior(SNP|Sj,Cov)j = Pr(Sj|SNP,Cov)*Prior(SNP|Cov) / ( Pr(Sj|Error,Cov)*Prior(Error|Cov) + Pr(Sj|SNP, Cov)*Prior(SNP|Cov) )
    Posterior(SNP|Sj,Cov)j: Posterior probability of a substitution site j to be a true SNP, when signal is Sj and variant at specific variant read coverage,cov.
    Pr(Sj|SNP,cov): Density distribution of Sj for SNPs
    Pr(Sj|Error): Density distribution of Sj for Errors
    Prior(SNP|Cov): Estimated SNP rate based variant read coverage
    Prior(Error|Cov): Estimated error rate based variant read coverage.
    Sj: Prior probability of a site j to be a true SNP.
  
Command line arguments
  Usage: ruby atlasSNPEvaluate_V3.0.rb -i [Substitution list] -x [option for XLR data] -e [priorerror] -o [output]
  -i is a required argument, which must be followed by the substitution list.
  -x is an optional argument. This select must be appended if the input is XLR data.
  -e is an optional argument with an associate value, Prior(Error|Cov) when variant coverage greater than 2. The default value is 0.1 if this option is ignored.
  -o is a required argument with an associate output file name with path.

Example:
  1. Evaluate on the XLR substitution list 'chr1.xm.SNP.list', with Prior(Error|Cov), 0.4 when variant coverage > 2. The
     output will be in file chr1/chr1.xm.SNP.list.eva.prior_0.4. 
       ruby atlasSNPEvaluate_V3.0.rb -i chr1.xm.SNP.list -x -e 0.4 -o chr1/chr1.xm.SNP.list.eva.prior_0.4
  2. Evaluate on the FLX substitution list 'chr1.xm.SNP.list', with default Prior(Error|Cov), 0.1 when variant coverage > 2.
     The output will be in file chr1/chr1.xm.SNP.list.eva.prior_0.1.
       ruby atlasSNPEvaluate_op_V3.1.rb -i chr1.xm.SNP.list -o chr1/chr1.xm.SNP.list.eva.prior_0.1"
  exit(2);
end

#Initialize the object of class 'AtlasEvaluate'
def initialize(optHash)
    @optHash=optHash
    setParameters()
end

#Logistic regression coefficients
$logistic_FLX={:intercept=>-5.36, :quality_score=>0.11, :swap=>-4.88, :NQS_pass=>1.07, :relpos=>0.91}
$logistic_XLR={:intercept=>-3.31, :quality_score=>0.1, :NQS_pass=>0.26, :swap=>-3.5, :relpos=>-0.37}

#Sj (variant read coverage <= 2)initialization
$errorPrior_test_less_3=[0.867, 0.867, 0.062, 0.062, 0.063, 0.063, 0.008, 0.008, 0.0003, 0.0003]
$snpPrior_test_less_3=[0.78, 0.78, 0.122, 0.122, 0.049, 0.049, 0.049, 0.049, 0.00001, 0.00001]

#Sj (variant read coverage >= 3)initialization
$errorPrior_test_greater_3=[0.676, 0.676, 0.102, 0.102, 0.074, 0.074, 0.028, 0.028, 0.12, 0.12]
$snpPrior_test_greater_3=[0.00001, 0.00001, 0.006, 0.006, 0.006, 0.006, 0.041, 0.041, 0.947, 0.947]

#Initialize the parameters
def setParameters()

if (!@optHash.key?('--input') or !@optHash.key?('--output'))
  AtlasEvaluate.usage("Option missing!")
  exit(2);
end

if @optHash.key?("--XLR")
  @XLR=true
  @logistic=$logistic_XLR
else
  @XLR=false
  @logistic=$logistic_FLX
end

if @optHash.key?("--priorerror")
  @perror_highcov=@optHash["--priorerror"]
  @psnp_highcov=1-@perror_highcov.to_f
else
  @perror_highcov=0.1
  @psnp_highcov=0.9
end

$stderr.puts "Prior(error|c): "+@perror_highcov.to_s
$stderr.puts "Prior(SNP|c):"+@psnp_highcov.to_s

#Create output file writer
@outputFile=@optHash["--output"]
@missedFile=@optHash["--output"]+".missed_reads"
@subsFile=@optHash["--input"]
end

#Evaluate the substitutions
def evaluate()
  @errPrior=[] 
  @snpPrior=[]
  perror=-10
  psnp=-10
  snp_cov=0
  outputWriter=File.new(@outputFile,"w")
  subsReader=File.open(@subsFile,"r")
  missedWriter=File.new(@missedFile,"w")
  outputWriter.print "refName\tcoordinate\trefBase\thomopolymer\trefEnv\tcoverage\tSNPBase\tadjustedQual\toriQual\tnumVariantReads\tnumAlternativeReads\treads_info\tarrayOfPr(error)i\tPrior(error|cov)j\tPrior(SNP|cov)j\tPr(Sj|error,cov)\tPr(Sj|SNP,cov)\tPrior(error|cov)\tPrior(SNP|cov)\tPosterior(error)\tPosterior(SNP)\tPosterior(SNP|Sj,Cov)j\n"
  subsReader.each do  |line|
  cols=line.split(/\s+/)

  pos, homo,cov,snpbase,num,readinfo =cols[1].to_i, cols[3].to_i,cols[5].to_i, cols[6], cols[9].to_i, cols[11]
  
  prior_error_a=1
  
  next if snpbase=="SNPBase"
  snp_cov=0
  
  lp_str=""
  readinfo.split(";").each do |r|
    if r=~ /^(\S+)\((\d+)\)(\S+)\((\d+)\)\S+([+|-])\S+\((\S+)\/(\S+)\/(\d+)\)\((\d+)\/(\d+)\)(\S+)/
      
      base=$1
        qual=$2.to_i
        name=$3
        dist=$4.to_i
        dir=$5
        sub=$6.to_f
        indel=$7.to_f
        tail=$8.to_i
        nqs=$9.to_i
        len=$10
        swap=$11
        relpos=dist/len.to_f
        
        next if base != snpbase
        next if nqs==nil or len==nil
        
        if base == 'N'
          logOdd = 0
          errorPosteriorFull = 1
        else
          if swap =="swap"
            swap = 1
            elsif swap =="mnp"
            swap = 0.5
            else
            swap = 0
          end
          
          #Logistic regression p-value
          predict = @logistic[:intercept] + @logistic[:quality_score]*qual + @logistic[:swap]*swap + @logistic[:NQS_pass]*nqs + @logistic[:relpos] * relpos
          logisticP = Math.exp(predict) / (1 + Math.exp(predict))
          
          pr_error_i=1-logisticP
          prior_error_a*=pr_error_i
          
        end
        
        snp_cov+=1
        
        lp_str<<pr_error_i.to_s+"|"
    end
  end
  
  if snp_cov==0
    missedWriter.puts line
  end
  
  #Select Sj
  if snp_cov>0
    outputWriter.print cols[0..-1].join("\t") + "\t" #Print the existing columns
    
    #Prior(SNP)j="Sj"
    prior_snp_a=1-BigDecimal(prior_error_a.to_s)
    bin = (BigDecimal(prior_snp_a.to_s) *10).floor #Determine the bin of Sa
  
    if bin==10
      bin=9
    end
  
    #Determine Sj sets according to variant reads coverage.
    if snp_cov<3
      @errPrior=$errorPrior_test_less_3 #Select Sj sets for Errors
      @snpPrior=$snpPrior_test_less_3 #Select Sj sets for SNPs
      perror=0.9 #Prior(Error|Cov)
      psnp=0.1 #Prior(SNP|Cov)
    elsif snp_cov>=3
      @errPrior=$errorPrior_test_greater_3
      @snpPrior=$snpPrior_test_greater_3
      
      perror=@perror_highcov.to_f
      psnp=@psnp_highcov.to_f
  end
  
  #Get Pr(Sj|Error,Cov) and Pr(Sj|SNP,Cov)
  ePrior = @errPrior[bin]
  sPrior = @snpPrior[bin]
  
  #For debugging
  if ePrior==nil or sPrior==nil
    $stderr.puts bin
    $stderr.puts prior_error_a
    $stderr.puts prior_snp_a
    $stderr.puts line
  end

  
  errorPosterior = ePrior * perror
  snpPosterior = sPrior * psnp
  
  #Get Posterior(SNP|Sj,Cov)j
  snpPosteriorFull = 1.0 / (1 + errorPosterior/snpPosterior)
 
  
  #print the useful variables and final Posterior(SNP|Sj,Cov)j
  outputWriter.print "#{lp_str}\t#{prior_error_a}\t#{prior_snp_a}\t#{ePrior}\t#{sPrior}\t#{perror}\t#{psnp}\t#{errorPosterior}\t#{snpPosterior}\t#{snpPosteriorFull}\n"
  end
  end
outputWriter.close
missedWriter.close
subsReader.close
end

end

optHash=AtlasEvaluate.processArguments()
evaluater=AtlasEvaluate.new(optHash)
evaluater.evaluate()

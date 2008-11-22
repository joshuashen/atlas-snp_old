$dir = ARGV[0]
$ref = ARGV[1]
$maq = File.expand_path(`which maq`.chomp)

def main
  dir = $dir
  ref = $ref
  flag = 1
  log = "Maq.log"
  maxMismatch = 3
  maxDist = 250
  pairedEndMapping(dir, ref, flag, log, maxMismatch, maxDist)
end

def pairedEndMapping(dir, ref, flag, log, maxMismatch, maxDist)
  # find pairs of bfq files, do maq
  bfqsFreads = Dir.entries(dir).sort.select {|file| file.match("bfq$") and file.match('F_reads')}
  
  bfqsFreads.each do |fReads|
    rReads = fReads.sub("F_reads", "R_reads")
    r1 = dir + '/' + fReads
    r2 = dir + '/' + rReads
    errorLog = r1 + '.stderr_log'
    stdoutlog = r1 + '.stdout_log'
    mapresult = r1 + '.map'
    maqstring = ''
    if File.exist?(r2)
      maqstring = "#{$maq} map  -n #{maxMismatch} -a #{maxDist}  #{mapresult}  #{ref} #{r1} #{r2} 2 >> #{log}"
      $stderr.puts "do maq: #{maqstring}"
    else
      $stderr.puts "Warning: cannot find the mate pair of #{fReads} :  #{rReads}"
    end
    
    if flag == 0 ## do locally
      $stderr.puts maqstring
      system(maqstring)
    elsif flag == 1 ## qsub
      shstring = '#!/bin/sh' + "\n" + '#$ -cwd' + "\n" + maqstring

      shf = File.new(r1+".sh", "w")
      shf.puts shstring
      shf.close
      system("qsub #{r1}.sh -e #{errorLog} -o #{stdoutlog}")
      
    elsif flag == 2 # bsub
      system("bsub -o #{stdoutlog} -e #{errorLog} #{maqstring}")
    end
  end
  return
end

main()    # do

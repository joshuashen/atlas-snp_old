# distribute the reads in a top-down fashion

fa = ARGV[0]

blatB = ARGV[1]

File.new(blatB, 'r').each do |line|
  

## 

str =  []
n = 0

while line = ARGF.gets do
  if line=~ /^(\d+)\s+/
    s = $1.to_i
    if s > n  # choose one
      puts str.sort_by {rand}[0]
      n = s
      str = []
    end
    str << line.chomp
  end
end


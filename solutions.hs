-- 'count' counts the occurrences of a given letter in a string
-- we build a counting function for each of letters A C G T, and run them on the input
dna :: String -> [Int]
dna xs = [f xs | f <- map count "ACGT"] where count c = length . filter (==c)

-- for each letter, if it's a T we change to U, else we leave it as is
rna :: String -> String
rna = map (\c -> if c=='T' then 'U' else c) . filter (`elem` "ACGT")

-- Reverse and complement, ignoring newlines in the input.
revc :: String -> String
revc = reverse . map compl . filter (/='\n')
  where compl 'A' = 'T'
        compl 'C' = 'G'
        compl 'G' = 'C'
        compl 'T' = 'A'
        compl _   = error "can only complement [ACGT]"

-- The probabilty of dominant phenotypes (one minus the probability of recessive homozygotes)
iprb :: [Double] -> Double
iprb [doms,hets,recs] = 1 - (0.25 * pr_hets*pr_hets + pr_hets*pr_recs + pr_recs*pr_recs)
  where total = doms + hets + recs
        pr_recs = recs/total
        pr_hets = hets/total

-- Prob of dominant phenotypes, from populations of pairs with given genotypes
iev :: [Double] -> Double
iev xs = let
  [pAA_AA,pAA_Aa,pAA_aa,pAa_Aa,pAa_aa,paa_aa] = map (/sum xs) xs
  in (pAA_AA + pAA_Aa + pAA_aa + 0.75*pAa_Aa + 0.5*pAa_aa)*2*sum xs

-- Hamming distance
hamm :: String -> String -> Int
hamm xs ys = length $ filter id $ zipWith (/=) xs ys

-- Find occrurrences of a string in another string
subs :: String -> String -> [Int]
subs str pat = go 1 str -- 1-indexed string, biologist was here
  where go i [] = []
        go i xs = let recurse = go (i+1) (tail xs) in if pat `isPrefixOf` xs then i : recurse else recurse

-- Number of rabbits after n iterations, each pair generates k new ones, and children are infertile for one iteration
fib :: Int -> Int -> Int
fib n k = let s = (0:1:zipWith (+) (map (k*) s) (tail s)) in s!!n

-- consensus sequence by majority rule, also output the frequency counts 
cons :: [String] -> [String]
cons xs = let ft = ftable xs in map consensus ft : showtable ft
   where
     consensus cs = snd . last . sort $ zip cs "ACGT"
     ftable xs | all null xs = []
               | otherwise = dna (map head xs) : ftable (map tail xs)
     showtable cols | all null cols = []
                    | otherwise = unwords (map (show . head) cols) : showtable (map tail cols)

-- calculate the log10 probability of a given string being generated, given GC-content
prob :: String -> [Double] -> [Double]
prob s = map (logBase 10 . prob1 s)
  where prob1 (x:xs) p = (if x `elem` "CG" then p/2 else (1-p)/2) * prob1 xs p
        prob1 [] _ = 1

-- convert mRNA to protein sequence
prot :: String -> String
prot = map lu . triples
  where lu x = case lookup x trans_tbl of
                Just y -> y
                Nothing -> error ("Unknown triple: "++show x)

-- helper function
triples :: String -> [String]
triples (a:b:c:rest) = [a,b,c]:triples rest
triples _ = []

-- find ORFs (also in revcompl, maybe that wasn't intented?)
orf :: String -> [String]
orf xs = let rs = revc xs
             candidates = map (prot . rna) [xs, tail xs, tail (tail xs), rs, tail rs, tail (tail rs)]
             orf1 ('M':rest) = ('M':takeWhile (/='!') rest) : orf1 rest
             orf1 (_:rest) = orf1 rest
             orf1 [] = []
         in concatMap orf1 candidates

-- remove introns from a sequence
splc :: [String] -> String
splc (s:introns) = prot . rna $ go s
  where go as | null as = []
              | otherwise = let matches = [ x | x <- introns, x `isPrefixOf` as ]
                            in case matches of [] -> head as : go (tail as)
                                               (x:_) -> go (drop (length x) as)

-- calculate the probability of at least one occurrence (=1-Pr(no occurrences))
-- of string 'str' in 'n' random strings with 'p' GC content
rstr :: Int -> Double -> String -> Double
rstr n p str = let [ps] = prob str [p]
               in 1-(1-10**ps)^n

-- Run them all on given inputs
main :: IO ()
main = do
  putStr "DNA: "
  putStrLn . unwords. map show . dna =<< readFile "bib2016/data/rosalind_dna.txt"
  putStr "\nRNA: "
  putStrLn . rna =<< readFile "bib2016/data/rosalind_rna.txt"
  putStr "\nREVC: "
  putStrLn . revc =<< readFile "bib2016/data/rosalind_revc.txt"
  putStr "\nIPRB: "
  putStrLn . show . iprb . map read . words =<< readFile "bib2016/data/rosalind_iprb.txt"
  putStr "\nIEV: "
  putStrLn . show . iev . map read . words =<< readFile "bib2016/data/rosalind_iev.txt"
  putStr "\nHAMM: "
  [l1,l2] <- lines `fmap` readFile "bib2016/data/rosalind_hamm.txt"
  putStrLn . show $ hamm l1 l2
  putStr "\nSUBS: "
  [l1,l2] <- lines `fmap` readFile "bib2016/data/rosalind_subs.txt"
  putStrLn $ unwords $ map show $ subs l1 l2
  putStr "\nFIB: "
  [n,k] <- map read `fmap` words `fmap` readFile "bib2016/data/rosalind_fib.txt"
  putStrLn . show $ fib n k
  putStrLn "\nCONS:"
  inp <- parseFasta "bib2016/data/rosalind_cons.txt"
  putStrLn . unlines $ cons inp
  [s,ps] <- lines `fmap` readFile "bib2016/data/rosalind_prob.txt"
  putStrLn . unwords . map show . prob s $ map read (words ps)
  putStr "\nPROT: "
  putStrLn . prot =<< readFile "bib2016/data/rosalind_prot.txt"
  putStrLn "\nORF: "
  [inp] <- parseFasta "bib2016/data/rosalind_orf.txt"
  putStrLn . unlines $ orf inp
  putStrLn "\nSPLC:"
  xs <- parseFasta "bib2016/data/rosalind_splc.txt"
  putStrLn $ splc xs
  putStrLn "\nRSTR:"
  [ns,ps,ss] <- words `fmap` readFile "bib2016/data/rosalind_rstr.txt"
  putStrLn $ show $ rstr (read ns) (read ps) ss

-- helper functions
  
trans_tbl :: [(String,Char)]
trans_tbl = [("UUU",'F'),("CUU",'L'),("AUU",'I'),("GUU",'V'),("UUC",'F'),("CUC",'L'),("AUC",'I'),("GUC",'V'),("UUA",'L'),("CUA",'L'),("AUA",'I'),("GUA",'V'),("UUG",'L'),("CUG",'L'),("AUG",'M'),("GUG",'V'),("UCU",'S'),("CCU",'P'),("ACU",'T'),("GCU",'A'),("UCC",'S'),("CCC",'P'),("ACC",'T'),("GCC",'A'),("UCA",'S'),("CCA",'P'),("ACA",'T'),("GCA",'A'),("UCG",'S'),("CCG",'P'),("ACG",'T'),("GCG",'A'),("UAU",'Y'),("CAU",'H'),("AAU",'N'),("GAU",'D'),("UAC",'Y'),("CAC",'H'),("AAC",'N'),("GAC",'D'),("UAA",'!'),("CAA",'Q'),("AAA",'K'),("GAA",'E'),("UAG",'!'),("CAG",'Q'),("AAG",'K'),("GAG",'E'),("UGU",'C'),("CGU",'R'),("AGU",'S'),("GGU",'G'),("UGC",'C'),("CGC",'R'),("AGC",'S'),("GGC",'G'),("UGA",'!'),("CGA",'R'),("AGA",'R'),("GGA",'G'),("UGG",'W'),("CGG",'R'),("AGG",'R'),("GGG",'G')]

-- reimplementations to avoid std library use (which is for wimps)
isPrefixOf [] _ = True
isPrefixOf _ [] = False
isPrefixOf (x:xs) (y:ys) = x==y && isPrefixOf xs ys

sort [] = []
sort (x:xs) = sort lt ++ [x] ++ sort gt
  where lt = filter (<= x) xs
        gt = filter (>x) xs
        
-- helper to read FASTA-formatted sequences
parseFasta :: FilePath -> IO [String]
parseFasta f = do
  inp <- lines `fmap` readFile f
  return (filter (not . null) $ go [] inp)
  where go acc (('>':_):ls) = concat (reverse acc) : go [] ls
        go acc (l:ls) = go (l:acc) ls
        go acc [] = [concat (reverse acc)]

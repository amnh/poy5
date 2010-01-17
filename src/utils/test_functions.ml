let tail = List.tl
let head = List.hd
let prepend a b = a :: b
let empty = []

let rec invert to_prepend pos len seq =
    if pos = 0 then
        if len = 0 then
            if to_prepend = empty then seq
            else 
                prepend (head to_prepend) (invert (tail to_prepend) 0 0 seq)
        else invert (prepend (head seq) to_prepend) 0 (pred len) (tail seq)
    else prepend (head seq) (invert empty (pred pos) len (tail seq))


let invert continuation pos len seq = 
    continuation (invert pos len seq)

let rec signed_inversion to_prepend pos len seq =
    if pos = 0 then
        if len = 0 then
            if to_prepend = empty then seq
            else prepend (not (head to_prepend)) (signed_inversion (tail to_prepend) 0 0
            seq)
        else signed_inversion (prepend (head seq) to_prepend) 0 (pred len) (tail seq)
    else prepend (head seq) (signed_inversion empty (pred pos) len (tail seq))

let rec cut chromosomes chrom_number position =
    if chrom_number = 0 then 
        if position = 0 then prepend empty chromosomes
        else
            let h = head (head chromosomes) in
            let chromosomes = 
                prepend (tail (head chromosomes)) (tail chromosomes) 
            in
            let chromosomes = cut chromosomes 0 (pred position) in
            prepend (prepend h  (head chromosomes)) (tail chromosomes)
    else prepend (head chromosomes) (cut (tail chromosomes) (pred chrom_number)
    position)

let rec full_invert res chrom = 
    if chrom = empty then res
    else full_invert (prepend (head chrom) res) (tail chrom)

let rec join chrom1 begin1 chrom2 begin2 =
    if chrom1 = empty then chrom2
    else if chrom2 = empty then chrom1 
    else
        if begin2 then
            if begin1 then 
                join (tail chrom1) begin1 (prepend (head chrom1) chrom2) begin2
            else prepend (head chrom1) (join (tail chrom1) begin1 chrom2 begin2)
        else join chrom1 begin1 (full_invert empty chrom2) true


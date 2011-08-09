
open Printf

let get_random_seq size = 
     let res = Array.make size 0 in
    for i = 0 to size-1 do
       res.(i) <- Random.int 4
    done;
    Array.iter (printf "%d,") res;
    let resstr = String.create size in
    let i = ref 0 in
    Array.iter (fun  x ->
        let _ = match x with
        | 0 ->  String.set resstr !i 'a'
        | 1 ->  String.set resstr !i 't'
        | 2 ->  String.set resstr !i 'c'
        | 3 ->  String.set resstr !i 'g'
        | _ ->  assert(false)
        in
        i := !i + 1;
    ) res;
    printf "\n{%s}\n%!" resstr;
    resstr

let ()  =
    Random.self_init();
    (*let line = input_line stdin in
    let size = (int_of_string line) in*)
    let size1 = Random.int 100 in
    let size1 = size1 + 1 in 
    printf "create random arr of size1 %d\n%!" size1;
    let resstr1 = get_random_seq size1 in 
    let size2 = Random.int 100 in
    let size2 = size2 + 1 in 
    printf "create random arr of size2 %d\n%!" size2;
    let resstr2 = get_random_seq size2 in 
    let fname = "ukkfile.in" in
    let oc =  open_out fname in
    fprintf oc "%s" resstr1;
    fprintf oc "\n";
    fprintf oc "%s" resstr2;
    fprintf oc "\n";
    close_out oc
    

open Printf

let debug_to_file = false
let debug_memsize = true

type dir = [`Start|`Insert|`Delete|`Align|`None] 
    
let print_intlst inlst oc =
    fprintf oc "[%!";
    List.iter (fprintf oc "%d,%!") inlst;
    fprintf oc "]\n%!"

let ukkonen_align seqarr1 seqarr2 cost_matrix =
    let debug = false and debug2 = false in
    let oc = 
        if debug_to_file then open_out "ukkss.out" 
        else stdout in
    let gapcode = Cost_matrix.Two_D.gap cost_matrix in
    (*why we don't have the min cost stored in Cost_matrix.Two_D? *)
    let pcm = Cost_matrix.Two_D.get_pure_cost_mat cost_matrix in
    let delta = 
        Array.fold_left (fun min arr ->
        let thismin = Array.fold_left (fun min cost ->
            if cost>0 && cost<min then cost else min
        ) (max_int/2) arr in
        if thismin<min then thismin else min
    ) (max_int/2) pcm
        (* let's get back to this later
    Array.fold_left (fun min arr ->
        let thismin = Array.fold_left (fun min cost ->
            if cost>0 && cost<min then cost else min
        ) (max_int/2) arr in
        if thismin<min then thismin else min
    ) (max_int/2) cost_matrix *)
    in
    let lenX = Array.length seqarr1 and lenY = Array.length seqarr2 in
    let seqarr1,seqarr2,lenX,lenY = 
        if lenX <= lenY then seqarr1,seqarr2,lenX,lenY 
        else seqarr2,seqarr1,lenY,lenX in
    let baseK = 0 in 
    let basebandwidth = lenY - lenX + 1 in
    (* (cost,possible_indel_number,dir)*)
    let base_ukkArr = Array.make (basebandwidth*lenX) (0,0,`None) in 
    (*we don't need this sizelst, it can be calculated, 
    * but it has maxlen=lenY and I'm lazzzzzzy*)
    let base_diagonal_sizelst =  (Array.to_list (Array.make basebandwidth lenX)) in
    if debug then begin
        fprintf oc "delta=%d,gapcode=%d\n%!" delta gapcode;
        if debug2 then fprintf oc "seqarr =\n%!";
    if debug2 then fprintf oc "[%!";
    if debug2 then Array.iter (fun x ->fprintf oc "%2i" x) seqarr1; 
    fprintf oc "]\n%!";
     if debug2 then fprintf oc "[%!";
    if debug2 then Array.iter (fun x ->fprintf oc "%2i" x) seqarr2; 
    if debug2 then fprintf oc "]\n%!";
    fprintf oc "ukkonen_align, lenX,lenY = %d,%d; baseK=%d,base band width =%d, \
    init ukkArr with size = %d\n%!" lenX lenY baseK basebandwidth (basebandwidth*lenX);
    if debug2 then print_intlst base_diagonal_sizelst oc ;
    end;
    (*get idx for ukkArr from normal (i,j),also return if (i.j) is on left/right border*)    
    let get_idx i j currentK diagonal_sizelst =
        let debug = false in
        assert(i>=0); assert(j>=0);
        let diagonal_before_me = j - (i-currentK) in
        let idx_in_my_diagonal = 
            if i<=j then i 
            else j 
        in
        let _,res = List.fold_left(fun (idx,acc) size ->
            if idx>diagonal_before_me then idx+1,acc
            else if idx=diagonal_before_me then idx+1,acc+idx_in_my_diagonal 
            else idx+1,acc+size
        ) (0,0) diagonal_sizelst in
        let at_leftborder =
            if diagonal_before_me = 0 then true
            else false
        and at_rightborder = 
            if diagonal_before_me = (List.length diagonal_sizelst)-1 then true
            else false
        in
        if debug then fprintf oc "[get_idx(%d,%d)->%d]\n%!" i j res;
        res,at_leftborder,at_rightborder
    in
    (*get cost from ukkArr*)
    let get_ukkcost idx ukkArr =
        let cost,max_gapnum,_ = ukkArr.(idx) in
        cost,max_gapnum
    in
    (*get cost between char of input seq*)
    let get_addcost i j =
        Cost_matrix.Two_D.cost seqarr1.(i) seqarr2.(j) cost_matrix
        (* let codei = seqarr1.(i) and codej = seqarr2.(j) in
        cost_matrix.(codei).(codej)*)
    in
    let get_gapcost insertGap i =
        if insertGap then
        Cost_matrix.Two_D.cost seqarr2.(i) gapcode cost_matrix
        else 
        Cost_matrix.Two_D.cost seqarr1.(i) gapcode cost_matrix
    in
    let print_ukkcell idx (cost,max_gapnum,dir) =
        fprintf oc "idx=%d:(cost=%d,max_g=%d," idx cost max_gapnum;
        match dir with
        | `Start -> fprintf oc "Start) \n%!"
        | `Insert -> fprintf oc "Insert) \n%!"
        | `Delete -> fprintf oc "Delete) \n%!"
        | `Align -> fprintf oc "Align) \n%!"
        | `None -> fprintf oc "None) \n%!"
    in
    let print_ukkArr inarr diagonal_sizelst k =
        fprintf oc "print ukkArr, k=%d" k;
        let i = ref (lenX-1) and j = ref (lenY-1) in
        while ((!i=0&& !j=0)=false) do
            let idx,_,_ = get_idx !i !j k diagonal_sizelst in
            let cost,max_gapnum,dir = inarr.(idx) in
            fprintf oc "(%d.%d):%!" !i !j;
            print_ukkcell idx (cost,max_gapnum,dir);
            match dir with
            | `Start -> 
                    assert(!i=0 && !j=0); 
            | `Insert -> 
                    j := !j - 1;
            | `Delete ->
                    i := !i - 1; 
            | `Align -> 
                    i := !i - 1; j := !j - 1;
            | `None -> fprintf oc "reach None cell!\n%!";
                    assert(false)
        done;
    in
    let ukk_traceback inarr diagonal_sizelst k =
        let k =
            if k>=lenX then lenX-1
            else k 
        in
        let debug = false in
        if debug then fprintf oc "ukk traceback, k=%d\n%!" k;
        let i = ref (lenX-1) and j = ref (lenY-1) in
        let reslst1 = ref [] and reslst2 = ref [] in
        while ((!i=0&& !j=0)=false) do
            let idx,_,_ = get_idx !i !j k diagonal_sizelst in
            let cost,max_gapnum,dir = inarr.(idx) in
            (*fprintf oc "(%d.%d):%!" !i !j;
            print_ukkcell idx (cost,max_gapnum,dir);*)
            match dir with
            | `Start -> 
                    assert(!i=0 && !j=0); 
            | `Insert -> 
                    reslst1 := gapcode :: !reslst1;
                    reslst2 := seqarr2.(!j) :: !reslst2;
                    j := !j - 1;
                    if !j<0 then 
                        print_ukkcell idx (cost,max_gapnum,dir);
            | `Delete ->
                    reslst1 := seqarr1.(!i) :: !reslst1;
                    reslst2 := gapcode :: !reslst2;
                    i := !i - 1; 
                    if !i<0 then
                        print_ukkcell idx (cost,max_gapnum,dir);
            | `Align ->
                    reslst1 := seqarr1.(!i) :: !reslst1;
                    reslst2 := seqarr2.(!j) :: !reslst2;
                    i := !i - 1; j := !j - 1;
                    if !i<0 || !j<0 then
                        print_ukkcell idx (cost,max_gapnum,dir);
            | `None -> failwith "reach None cell!"
        done;
        if debug then print_intlst !reslst1 oc;
        if debug then print_intlst !reslst2 oc;
        Array.of_list !reslst1, Array.of_list !reslst2
    in
   (*function for updating a cell in ukkArr start*)
(*update non-border cell, with i-1,j-1 and i-1,j and i,j-1*)
    let update_internal_cell ukkidx idxi idxj currentK diagonal_sizelst ukkArr =
        let debug =false in
        if debug then  fprintf oc "update internal cell ->%!";
        let costfromL,max_gapnum_fromL =
            if idxj=0 then max_int/2,0
            else
                let ukkidxL,_,_ = get_idx idxi (idxj-1) currentK diagonal_sizelst in
                get_ukkcost ukkidxL ukkArr
        in
        let costfromR,max_gapnum_fromR = 
            if idxi=0 then max_int/2,0
            else
            let ukkidxR,_,_ = get_idx (idxi-1) idxj currentK diagonal_sizelst in
            get_ukkcost ukkidxR ukkArr in
        let costfromM,max_gapnum_fromM = 
            if idxi=0||idxj=0 then max_int/2,0
            else
            let ukkidxM,_,_ = get_idx (idxi-1) (idxj-1) currentK
            diagonal_sizelst in
            get_ukkcost ukkidxM ukkArr in
        if debug then  fprintf oc "costfromL=%d,costfromR=%d,costfromM=%d---> res="
        costfromL costfromR costfromM ;
        let costL = costfromL + get_gapcost true idxj in
        let costR = costfromR + get_gapcost false idxi in
        let costM = costfromM + get_addcost idxi idxj in
        if costL<=costR && costL<=costM then
            ukkArr.(ukkidx) <- (costL,max_gapnum_fromL+1,`Insert)
        else if costR<=costL && costR<=costM then
            ukkArr.(ukkidx) <- (costR,max_gapnum_fromR+1,`Delete)
        else
            let max_gapnum_fromM = if costM=costfromM then max_gapnum_fromM else
                max_gapnum_fromM + 1 in
            ukkArr.(ukkidx) <- (costM,max_gapnum_fromM,`Align)
        ;
        if debug then print_ukkcell ukkidx ukkArr.(ukkidx);
    in
(*update left border cell.i.j, with i-1,j-1 and i-1,j*)
    let update_left_border_cell ukkidx idxi idxj currentK diagonal_sizelst ukkArr = 
        let debug =false in
        if debug then  fprintf oc "update left border cell -> %!";
        let costfromR,max_gapnum_fromR = 
            if idxi=0 then max_int/2,0
            else
            let ukkidxR,_,_ = get_idx (idxi-1) idxj currentK diagonal_sizelst in
            get_ukkcost ukkidxR ukkArr in
        let costfromM,max_gapnum_fromM = 
            if idxi=0||idxj=0 then max_int/2,0
            else
            let ukkidxM,at_leftborder,_ = get_idx (idxi-1) (idxj-1) currentK diagonal_sizelst in
            assert(at_leftborder);
            get_ukkcost ukkidxM ukkArr in
        if debug then  fprintf oc "costfromR=%d,costfromM=%d  ---> res=" costfromR costfromM;
        let costR = costfromR + get_gapcost false idxi in
        let costM = costfromM + get_addcost idxi idxj in
        if costR<=costM then
            ukkArr.(ukkidx) <- (costR,max_gapnum_fromR,`Delete)
        else
            let max_gapnum_fromM = if costM=costfromM then max_gapnum_fromM else
                max_gapnum_fromM + 1 in
            ukkArr.(ukkidx) <- (costM,max_gapnum_fromM,`Align)
        ;
        if debug then print_ukkcell ukkidx ukkArr.(ukkidx);
    in
(*update right border cell.i.j with i-1,j-1 and i,j-1*)
    let update_right_border_cell ukkidx idxi idxj currentK diagonal_sizelst ukkArr =
        let debug =false in
        if debug then  fprintf oc "update rigght border cell ->%!";
        let costfromL,max_gapnum_fromL = 
            if idxj=0 then max_int/2,0
            else
            let ukkidxL,_,_ = get_idx idxi (idxj-1) currentK diagonal_sizelst in
            get_ukkcost ukkidxL ukkArr in
        let costfromM,max_gapnum_fromM = 
            if idxi=0||idxj=0 then max_int/2,0
            else
            let ukkidxM,_,at_rightborder = get_idx (idxi-1) (idxj-1) currentK diagonal_sizelst in
            assert(at_rightborder);
            get_ukkcost ukkidxM ukkArr in
        if debug then  fprintf oc "costfromL=%d,costfromM=%d  ---> res=" costfromL costfromM; 
        let costL = costfromL + get_gapcost true idxj in
        let costM = costfromM + get_addcost idxi idxj in
        if costL<=costM then
            ukkArr.(ukkidx) <- (costL,max_gapnum_fromL+1,`Insert)
        else
            let max_gapnum_fromM = if costM=costfromM then max_gapnum_fromM else
                max_gapnum_fromM + 1 in
            ukkArr.(ukkidx) <- (costM,max_gapnum_fromM,`Align)
        ;
        if debug then print_ukkcell ukkidx ukkArr.(ukkidx);
    in
    let update_central_diagonal_cell ukkidx idxi idxj currentK diagonal_sizelst ukkArr =
        let debug = false in
        let costfromM,max_gapnum_fromM = 
            if idxi=0||idxj=0 then max_int/2,0
            else
            let ukkidxM,_,at_rightborder = get_idx (idxi-1) (idxj-1) currentK diagonal_sizelst in
            assert(at_rightborder);
            get_ukkcost ukkidxM ukkArr in
        let costM = costfromM + get_addcost idxi idxj in
        let max_gapnum_fromM = if costM=costfromM then max_gapnum_fromM else
                max_gapnum_fromM + 1 in
        ukkArr.(ukkidx) <- (costM,max_gapnum_fromM,`Align);
        if debug then print_ukkcell ukkidx ukkArr.(ukkidx);
    in
    let update_a_cell idxi idxj currentK diagonal_sizelst ukkArr =
        let debug = false in
        if debug then fprintf oc "\n ++++++++++++ update cell.%d.%d,currentK=%d -> \n%!" idxi idxj currentK;
        if idxi=0&&idxj=0 then ()(*ukkArr.(0) is set to 0 always.*)
        else begin
        let ukkidx,at_leftborder,at_rightborder = get_idx idxi idxj currentK diagonal_sizelst in
        if at_leftborder&&at_rightborder then (*this happens when lenX=lenY*)
            update_central_diagonal_cell ukkidx idxi idxj currentK
            diagonal_sizelst ukkArr
        else if (at_leftborder) then
             update_left_border_cell ukkidx idxi idxj currentK diagonal_sizelst
             ukkArr
        else if (at_rightborder) then
            update_right_border_cell ukkidx idxi idxj currentK diagonal_sizelst
            ukkArr
        else
            update_internal_cell ukkidx idxi idxj currentK diagonal_sizelst ukkArr
        ;
        end
    in
(*function for increaseing ukkArr size and diagonal_sizelst.
* initial case: when oldk=newk, nothing happen*)
    let update_diagonal_sizelst inlst oldk newk inukkarr =
        let addsize = ref(List.hd inlst) in
        let addlst = ref [] in
        let resukkarr = ref inukkarr in
        for diffk = 1 to newk-oldk do
            addsize := !addsize -1 ;
            addlst := !addsize :: !addlst;
            resukkarr :=
                Array_ops.array_append (Array.make !addsize (0,0,`None)) !resukkarr;
            resukkarr := 
                Array_ops.array_append !resukkarr (Array.make !addsize (0,0,`None));
        done;
        let addlst = !addlst in
        if debug_memsize then 
            Printf.printf "ukkArr size set to %d\n%!" (2*(List.length
            addlst)+(List.length inlst));
        addlst@inlst@(List.rev addlst),!resukkarr
    in
(*test if we reject the currentK*)
    let ukktest thresholdT p diagonal_sizelst ukkArr =
        let debug = false in
        let newk = 
            if p>=lenX then lenX-1
            else p in
        let oldk = 
            ((List.length diagonal_sizelst) - basebandwidth)/2 
        in 
        let new_diagonal_sizelst,new_ukkArr = update_diagonal_sizelst
        diagonal_sizelst oldk newk ukkArr in
        if debug then 
        fprintf oc "\n ukktest on \
        newT=%d,newp=%d,newk=%d(oldk=%d),diagonal_sizelst(size=%d)= %!"
        thresholdT p newk oldk (List.length new_diagonal_sizelst);
        if debug then print_intlst new_diagonal_sizelst oc;
        (*cell in rectangle area of [lefttop,leftbottom,rightbottom,righttop] will not change*)
        let in_non_change_zone i j =
            let lefttopi,lefttopj = 0,0 
            and leftbottomi,leftbottomj = oldk,0
            and righttopi,righttopj = 0,lenY-lenX+oldk in
            ( i>=lefttopi && i<= leftbottomi && j>= lefttopj && j<=righttopj)
        in
        let outside_diagonal_area i j = 
            ( (i>j)&&(i-j)>newk ) || ( (i<j)&&(j-i)>(lenY-lenX+newk) )
        in
        for i=0 to lenX-1 do
            for j= 0 to lenY-1 do
                if (in_non_change_zone i j)  || (outside_diagonal_area i j) 
                  (*  (i=0&&j=0) || ( (i>j)&&(i-j)>newk ) || ( (i<j)&&(j-i)>(lenY-lenX+newk) ) 
                  *)then ()
                else 
                update_a_cell i j newk new_diagonal_sizelst new_ukkArr;
            done;
        done;
        let ukkidx,_,_ = get_idx (lenX-1) (lenY-1) newk new_diagonal_sizelst in 
        let cost,max_gapnum,_ = new_ukkArr.(ukkidx) in
        if debug then 
            fprintf oc "update ukkArr done,check ukkArr.%d.%d\n%!" (lenX-1) (lenY-1);
        if debug then print_ukkcell ukkidx new_ukkArr.(ukkidx);
        if debug then print_ukkArr new_ukkArr new_diagonal_sizelst newk;
        (*according to ukkonen's paper, we test if current cost is smaller than
        * thresholdT, if so , we have our result, if not, double T and do it again*)
        cost<=thresholdT,
        cost,
        max_gapnum,
        new_diagonal_sizelst,
        new_ukkArr
    in
    let rec increaseT newT diagonal_sizelst ukkArr =
        let debug = false in
        let p = ( newT - (lenY-lenX))/2 in
        if debug then fprintf oc "increaseT, newT=%d,p=%d,%!" newT p;
        let testres,new_cost,new_max_gapnum,new_diagonal_sizelst,new_ukkArr= 
            ukktest newT p diagonal_sizelst ukkArr
        in
        if debug then fprintf oc "cost=%d\n%!" new_cost;
        let newp = (newT*2 -(lenY-lenX))/2 in
        if testres=false then 
            if (new_cost < newT*2)||(newp>new_max_gapnum) then begin
                if debug then fprintf oc "cost=%d,end by K(%d>%d) or T(%d<%d),ukk trace back\n%!" 
                new_cost newp new_max_gapnum new_cost (newT*2); 
                let alied_seqarr1,alied_seqarr2 = 
                    ukk_traceback new_ukkArr new_diagonal_sizelst p in
                if debug_to_file then close_out oc;
                alied_seqarr1,alied_seqarr2,new_cost
            end
            else
            increaseT (newT*2) new_diagonal_sizelst new_ukkArr
        else begin
            let alied_seqarr1,alied_seqarr2 = 
                ukk_traceback new_ukkArr new_diagonal_sizelst p in
            if debug then fprintf oc "cost=%d,end by T,ukk trace back\n%!" new_cost;
            if debug_to_file then close_out oc;
            alied_seqarr1,alied_seqarr2,new_cost
        end
        ;
    in
    (*init the first row of baseband*)
    base_ukkArr.(0) <- (0,0,`Start);
    for i = 1 to basebandwidth-1 do
        let debug = false in
        let preukkidx,_,_ = get_idx 0 (i-1) baseK base_diagonal_sizelst in
        let precost,_ = get_ukkcost preukkidx base_ukkArr in
        let thiscost = get_gapcost true i in
        let ukkidx,_,_ = get_idx 0 i baseK base_diagonal_sizelst in
        base_ukkArr.(ukkidx) <- (thiscost+precost,i,`Insert);
        if debug then print_ukkcell ukkidx base_ukkArr.(ukkidx);
    done;
    increaseT (basebandwidth*delta) base_diagonal_sizelst base_ukkArr
(*
let () =
    let oc =  open_out "ukkss.out" in
    fprintf oc "delta=%d\n%!" delta;
    let filename = "ukkfile.in" in
    let lines = ref [] in
    let chan = open_in filename in
    try 
    while true do
    lines := input_line chan :: !lines
    done;
    with End_of_file ->
    close_in chan;
    let lines = List.rev !lines in
    fprintf oc "read in str=\n%!";
    List.iter (fun str -> fprintf oc "%s\n%!" str) lines;
    let seqlst = List.fold_right (fun seqstr acc ->
        if (String.length seqstr)>0 && 
        (seqstr.[0]='a'||seqstr.[0]='c'||seqstr.[0]='t'||seqstr.[0]='g') then
        begin
        let seqlst = ref [] in
        for i=0 to (String.length seqstr)-1 do
            seqlst := 
            (match seqstr.[i] with
            | 'a' -> 0::!seqlst
            | 'c' -> 1::!seqlst
            | 'g' -> 2::!seqlst
            | 't' -> 3::!seqlst
            | _ -> !seqlst )
            ;
        done;
        (Array.of_list (List.rev !seqlst))::acc
        end
        else acc
    ) lines [[||]] 
    in
    fprintf oc "seqarr =\n%!";
     List.iter (
        fun seqlst ->
            fprintf oc "[%!";
            Array.iter (fun x ->fprintf oc "%d" x) seqlst; 
        fprintf oc "]\n%!";
        ) seqlst;
    let seqarr1,seqarr2 = List.hd seqlst, List.nth seqlst 1 in
    let seqarr1,seqarr2 = 
        Array.append [|0|] seqarr1, 
        Array.append [|0|] seqarr2
    in
    let seqarr1,seqarr2 = 
        if (Array.length seqarr1)<=(Array.length seqarr2) then
            seqarr1,seqarr2
        else 
            seqarr2,seqarr1
    in
    ukkonen_align seqarr1 seqarr2 oc;
    close_out oc
*)

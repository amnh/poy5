(* A List of unicode equivalents for some of the PDF ZapfDingbats Encoding.
There is no standard source for this. *)
let dingbatmap = [ "/a1", [0x2701]; "/a2", [0x2702]; "/a202", [0x2703]; "/a3",
[0x2704]; "/a4", [0x260e]; "/a5", [0x2706]; "/a119", [0x2707]; "/a118",
[0x2708]; "/a117", [0x2709]; "/a11", [0x261B]; "/a12", [0x261E]; "/a13",
[0x270C]; "/a14", [0x270D]; "/a15", [0x270E]; "/a16", [0x270F]; "/a105",
[0x2710]; "/a17", [0x2711]; "/a18", [0x2712]; "/a19", [0x2713]; "/a20",
[0x2714]; "/a21", [0x2715]; "/a22", [0x2716]; "/a23", [0x2717]; "/a24",
[0x2718]; "/a25", [0x2719]; "/a26", [0x271A]; "/a27", [0x271B]; "/a28",
[0x271C]; "/a6", [0x271D]; "/a7", [0x271E]; "/a8", [0x271F]; "/a9", [0x2720];
"/a10", [0x2721]; "/a29", [0x2722]; "/a30", [0x2723]; "/a31", [0x2724]; "/a32",
[0x2725]; "/a33", [0x2726]; "/a34", [0x2727]; "/a35", [0x2605]; "/a36",
[0x2606]; "/a37", [0x272A]; "/a38", [0x272B]; "/a39", [0x272C]; "/a40",
[0x272D]; "/a41", [0x272E]; "/a42", [0x272F]; "/a43", [0x2730]; "/a44",
[0x2731]; "/a45", [0x2732]; "/a46", [0x2733]; "/a47", [0x2734]; "/a48",
[0x2735]; "/a49", [0x2736]; "/a50", [0x2737]; "/a51", [0x2738]; "/a52",
[0x2739]; "/a53", [0x273A]; "/a54", [0x273B]; "/a55", [0x273C]; "/a56",
[0x273D]; "/a57", [0x273E]; "/a58", [0x273F]; "/a59", [0x2740]; "/a60",
[0x2741]; "/a61", [0x2742]; "/a62", [0x2743]; "/a63", [0x2744]; "/a64",
[0x2745]; "/a65", [0x2746]; "/a66", [0x2747]; "/a67", [0x2748]; "/a68",
[0x2749]; "/a69", [0x274A]; "/a70", [0x274B]; "/a71", [0x25CF]; "/a72",
[0x274D]; "/a73", [0x25A0]; "/a74", [0x274F]; "/a203", [0x2750]; "/a75",
[0x2751]; "/a204", [0x2752]; "/a76", [0x25B2]; "/a77", [0x25BC]; "/a78",
[0x25C6]; "/a79", [0x2756]; "/a81", [0x25D7]; "/a82", [0x2758]; "/a83",
[0x2759]; "/a84", [0x275A]; "/a97", [0x275B]; "/a98", [0x275C]; "/a99",
[0x275D]; "/a100", [0x275E]; "/a101", [0x2761]; "/a102", [0x2762]; "/a103",
[0x2763]; "/a104", [0x2764]; "/a106", [0x2765]; "/a107", [0x2766]; "/a108",
[0x2767]; "/a112", [0x2663]; "/a111", [0x2666]; "/a110", [0x2665]; "/a109",
[0x2660]; "/a120", [0x2780]; "/a121", [0x2781]; "/a122", [0x2782]; "/a123",
[0x2783]; "/a124", [0x2784]; "/a125", [0x2785]; "/a126", [0x2786]; "/a127",
[0x2787]; "/a128", [0x2788]; "/a129", [0x2789]; "/a130", [0x2776]; "/a131",
[0x2777]; "/a132", [0x2778]; "/a133", [0x2779]; "/a134", [0x277A]; "/a135",
[0x277B]; "/a136", [0x277C]; "/a137", [0x277D]; "/a138", [0x277E]; "/a139",
[0x277F]; "/a140", [0x2780]; "/a141", [0x2781]; "/a142", [0x2782]; "/a143",
[0x2783]; "/a144", [0x2784]; "/a145", [0x2785]; "/a146", [0x2786]; "/a147",
[0x2787]; "/a148", [0x2788]; "/a149", [0x2789]; "/a150", [0x278A]; "/a151",
[0x278B]; "/a152", [0x278C]; "/a153", [0x278D]; "/a154", [0x278E]; "/a155",
[0x278F]; "/a156", [0x2790]; "/a157", [0x2791]; "/a158", [0x2792]; "/a159",
[0x2793]; "/a160", [0x2794]; "/a161", [0x279D];
"/a163", [0x2194]; (* Not an exact match *)
"/a164", [0x2195]; (* Ditto *)
"/a196", [0x2798]; "/a165", [0x2799]; "/a192", [0x279A]; "/a166", [0x279B];
"/a167", [0x279C]; "/a168", [0x279D]; "/a169", [0x279E]; "/a170", [0x279F];
"/a171", [0x27A0]; "/a172", [0x27A1]; "/a173", [0x27A2]; "/a162", [0x27A3];
"/a174", [0x27A4]; "/a175", [0x27A5]; "/a176", [0x27A6]; "/a177", [0x27A7];
"/a178", [0x27A8]; "/a179", [0x27A9]; "/a193", [0x27AA]; "/a180", [0x27AB];
"/a199", [0x27AC]; "/a181", [0x27AD]; "/a200", [0x27AE]; "/a182", [0x27AF];
"/a201", [0x2781]; "/a183", [0x27B2]; "/a184", [0x27B3]; "/a197", [0x27B4];
"/a185", [0x27B5]; "/a194", [0x27B6]; "/a198", [0x27B7]; "/a186", [0x27B8];
"/a195", [0x27B9]; "/a187", [0x27BA]; "/a188", [0x27BB]; "/a189", [0x27BC];
"/a190", [0x27BD]; "/a191", [0x27BE]]

(* The Adobe Glyph list, 2.0 *)
let glyphmap = ["/A", [0x0041]; "/AE", [0x00C6]; "/AEacute", [0x01FC];
"/AEmacron", [0x01E2]; "/AEsmall", [0xF7E6]; "/Aacute", [0x00C1];
"/Aacutesmall", [0xF7E1]; "/Abreve", [0x0102]; "/Abreveacute", [0x1EAE];
"/Abrevecyrillic", [0x04D0]; "/Abrevedotbelow", [0x1EB6]; "/Abrevegrave",
[0x1EB0]; "/Abrevehookabove", [0x1EB2]; "/Abrevetilde", [0x1EB4]; "/Acaron",
[0x01CD]; "/Acircle", [0x24B6]; "/Acircumflex", [0x00C2]; "/Acircumflexacute",
[0x1EA4]; "/Acircumflexdotbelow", [0x1EAC]; "/Acircumflexgrave", [0x1EA6];
"/Acircumflexhookabove", [0x1EA8]; "/Acircumflexsmall", [0xF7E2];
"/Acircumflextilde", [0x1EAA]; "/Acute", [0xF6C9]; "/Acutesmall", [0xF7B4];
"/Acyrillic", [0x0410]; "/Adblgrave", [0x0200]; "/Adieresis", [0x00C4];
"/Adieresiscyrillic", [0x04D2]; "/Adieresismacron", [0x01DE]; "/Adieresissmall",
[0xF7E4]; "/Adotbelow", [0x1EA0]; "/Adotmacron", [0x01E0]; "/Agrave", [0x00C0];
"/Agravesmall", [0xF7E0]; "/Ahookabove", [0x1EA2]; "/Aiecyrillic", [0x04D4];
"/Ainvertedbreve", [0x0202]; "/Alpha", [0x0391]; "/Alphatonos", [0x0386];
"/Amacron", [0x0100]; "/Amonospace", [0xFF21]; "/Aogonek", [0x0104]; "/Aring",
[0x00C5]; "/Aringacute", [0x01FA]; "/Aringbelow", [0x1E00]; "/Aringsmall",
[0xF7E5]; "/Asmall", [0xF761]; "/Atilde", [0x00C3]; "/Atildesmall", [0xF7E3];
"/Aybarmenian", [0x0531]; "/B", [0x0042]; "/Bcircle", [0x24B7]; "/Bdotaccent",
[0x1E02]; "/Bdotbelow", [0x1E04]; "/Becyrillic", [0x0411]; "/Benarmenian",
[0x0532]; "/Beta", [0x0392]; "/Bhook", [0x0181]; "/Blinebelow", [0x1E06];
"/Bmonospace", [0xFF22]; "/Brevesmall", [0xF6F4]; "/Bsmall", [0xF762];
"/Btopbar", [0x0182]; "/C", [0x0043]; "/Caarmenian", [0x053E]; "/Cacute",
[0x0106]; "/Caron", [0xF6CA]; "/Caronsmall", [0xF6F5]; "/Ccaron", [0x010C];
"/Ccedilla", [0x00C7]; "/Ccedillaacute", [0x1E08]; "/Ccedillasmall", [0xF7E7];
"/Ccircle", [0x24B8]; "/Ccircumflex", [0x0108]; "/Cdot", [0x010A];
"/Cdotaccent", [0x010A]; "/Cedillasmall", [0xF7B8]; "/Chaarmenian", [0x0549];
"/Cheabkhasiancyrillic", [0x04BC]; "/Checyrillic", [0x0427];
"/Chedescenderabkhasiancyrillic", [0x04BE]; "/Chedescendercyrillic", [0x04B6];
"/Chedieresiscyrillic", [0x04F4]; "/Cheharmenian", [0x0543];
"/Chekhakassiancyrillic", [0x04CB]; "/Cheverticalstrokecyrillic", [0x04B8];
"/Chi", [0x03A7]; "/Chook", [0x0187]; "/Circumflexsmall", [0xF6F6];
"/Cmonospace", [0xFF23]; "/Coarmenian", [0x0551]; "/Csmall", [0xF763]; "/D",
[0x0044]; "/DZ", [0x01F1]; "/DZcaron", [0x01C4]; "/Daarmenian", [0x0534];
"/Dafrican", [0x0189]; "/Dcaron", [0x010E]; "/Dcedilla", [0x1E10]; "/Dcircle",
[0x24B9]; "/Dcircumflexbelow", [0x1E12]; "/Dcroat", [0x0110]; "/Ddotaccent",
[0x1E0A]; "/Ddotbelow", [0x1E0C]; "/Decyrillic", [0x0414]; "/Deicoptic",
[0x03EE]; "/Delta", [0x2206]; "/Deltagreek", [0x0394]; "/Dhook", [0x018A];
"/Dieresis", [0xF6CB]; "/DieresisAcute", [0xF6CC]; "/DieresisGrave", [0xF6CD];
"/Dieresissmall", [0xF7A8]; "/Digammagreek", [0x03DC]; "/Djecyrillic", [0x0402];
"/Dlinebelow", [0x1E0E]; "/Dmonospace", [0xFF24]; "/Dotaccentsmall", [0xF6F7];
"/Dslash", [0x0110]; "/Dsmall", [0xF764]; "/Dtopbar", [0x018B]; "/Dz", [0x01F2];
"/Dzcaron", [0x01C5]; "/Dzeabkhasiancyrillic", [0x04E0]; "/Dzecyrillic",
[0x0405]; "/Dzhecyrillic", [0x040F]; "/E", [0x0045]; "/Eacute", [0x00C9];
"/Eacutesmall", [0xF7E9]; "/Ebreve", [0x0114]; "/Ecaron", [0x011A];
"/Ecedillabreve", [0x1E1C]; "/Echarmenian", [0x0535]; "/Ecircle", [0x24BA];
"/Ecircumflex", [0x00CA]; "/Ecircumflexacute", [0x1EBE]; "/Ecircumflexbelow",
[0x1E18]; "/Ecircumflexdotbelow", [0x1EC6]; "/Ecircumflexgrave", [0x1EC0];
"/Ecircumflexhookabove", [0x1EC2]; "/Ecircumflexsmall", [0xF7EA];
"/Ecircumflextilde", [0x1EC4]; "/Ecyrillic", [0x0404]; "/Edblgrave", [0x0204];
"/Edieresis", [0x00CB]; "/Edieresissmall", [0xF7EB]; "/Edot", [0x0116];
"/Edotaccent", [0x0116]; "/Edotbelow", [0x1EB8]; "/Efcyrillic", [0x0424];
"/Egrave", [0x00C8]; "/Egravesmall", [0xF7E8]; "/Eharmenian", [0x0537];
"/Ehookabove", [0x1EBA]; "/Eightroman", [0x2167]; "/Einvertedbreve", [0x0206];
"/Eiotifiedcyrillic", [0x0464]; "/Elcyrillic", [0x041B]; "/Elevenroman",
[0x216A]; "/Emacron", [0x0112]; "/Emacronacute", [0x1E16]; "/Emacrongrave",
[0x1E14]; "/Emcyrillic", [0x041C]; "/Emonospace", [0xFF25]; "/Encyrillic",
[0x041D]; "/Endescendercyrillic", [0x04A2]; "/Eng", [0x014A]; "/Enghecyrillic",
[0x04A4]; "/Enhookcyrillic", [0x04C7]; "/Eogonek", [0x0118]; "/Eopen", [0x0190];
"/Epsilon", [0x0395]; "/Epsilontonos", [0x0388]; "/Ercyrillic", [0x0420];
"/Ereversed", [0x018E]; "/Ereversedcyrillic", [0x042D]; "/Escyrillic", [0x0421];
"/Esdescendercyrillic", [0x04AA]; "/Esh", [0x01A9]; "/Esmall", [0xF765]; "/Eta",
[0x0397]; "/Etarmenian", [0x0538]; "/Etatonos", [0x0389]; "/Eth", [0x00D0];
"/Ethsmall", [0xF7F0]; "/Etilde", [0x1EBC]; "/Etildebelow", [0x1E1A]; "/Euro",
[0x20AC]; "/Ezh", [0x01B7]; "/Ezhcaron", [0x01EE]; "/Ezhreversed", [0x01B8];
"/F", [0x0046]; "/Fcircle", [0x24BB]; "/Fdotaccent", [0x1E1E]; "/Feharmenian",
[0x0556]; "/Feicoptic", [0x03E4]; "/Fhook", [0x0191]; "/Fitacyrillic", [0x0472];
"/Fiveroman", [0x2164]; "/Fmonospace", [0xFF26]; "/Fourroman", [0x2163];
"/Fsmall", [0xF766]; "/G", [0x0047]; "/GBsquare", [0x3387]; "/Gacute", [0x01F4];
"/Gamma", [0x0393]; "/Gammaafrican", [0x0194]; "/Gangiacoptic", [0x03EA];
"/Gbreve", [0x011E]; "/Gcaron", [0x01E6]; "/Gcedilla", [0x0122]; "/Gcircle",
[0x24BC]; "/Gcircumflex", [0x011C]; "/Gcommaaccent", [0x0122]; "/Gdot",
[0x0120]; "/Gdotaccent", [0x0120]; "/Gecyrillic", [0x0413]; "/Ghadarmenian",
[0x0542]; "/Ghemiddlehookcyrillic", [0x0494]; "/Ghestrokecyrillic", [0x0492];
"/Gheupturncyrillic", [0x0490]; "/Ghook", [0x0193]; "/Gimarmenian", [0x0533];
"/Gjecyrillic", [0x0403]; "/Gmacron", [0x1E20]; "/Gmonospace", [0xFF27];
"/Grave", [0xF6CE]; "/Gravesmall", [0xF760]; "/Gsmall", [0xF767]; "/Gsmallhook",
[0x029B]; "/Gstroke", [0x01E4]; "/H", [0x0048]; "/H18533", [0x25CF]; "/H18543",
[0x25AA]; "/H18551", [0x25AB]; "/H22073", [0x25A1]; "/HPsquare", [0x33CB];
"/Haabkhasiancyrillic", [0x04A8]; "/Hadescendercyrillic", [0x04B2];
"/Hardsigncyrillic", [0x042A]; "/Hbar", [0x0126]; "/Hbrevebelow", [0x1E2A];
"/Hcedilla", [0x1E28]; "/Hcircle", [0x24BD]; "/Hcircumflex", [0x0124];
"/Hdieresis", [0x1E26]; "/Hdotaccent", [0x1E22]; "/Hdotbelow", [0x1E24];
"/Hmonospace", [0xFF28]; "/Hoarmenian", [0x0540]; "/Horicoptic", [0x03E8];
"/Hsmall", [0xF768]; "/Hungarumlaut", [0xF6CF]; "/Hungarumlautsmall", [0xF6F8];
"/Hzsquare", [0x3390]; "/I", [0x0049]; "/IAcyrillic", [0x042F]; "/IJ", [0x0132];
"/IUcyrillic", [0x042E]; "/Iacute", [0x00CD]; "/Iacutesmall", [0xF7ED];
"/Ibreve", [0x012C]; "/Icaron", [0x01CF]; "/Icircle", [0x24BE]; "/Icircumflex",
[0x00CE]; "/Icircumflexsmall", [0xF7EE]; "/Icyrillic", [0x0406]; "/Idblgrave",
[0x0208]; "/Idieresis", [0x00CF]; "/Idieresisacute", [0x1E2E];
"/Idieresiscyrillic", [0x04E4]; "/Idieresissmall", [0xF7EF]; "/Idot", [0x0130];
"/Idotaccent", [0x0130]; "/Idotbelow", [0x1ECA]; "/Iebrevecyrillic", [0x04D6];
"/Iecyrillic", [0x0415]; "/Ifraktur", [0x2111]; "/Igrave", [0x00CC];
"/Igravesmall", [0xF7EC]; "/Ihookabove", [0x1EC8]; "/Iicyrillic", [0x0418];
"/Iinvertedbreve", [0x020A]; "/Iishortcyrillic", [0x0419]; "/Imacron", [0x012A];
"/Imacroncyrillic", [0x04E2]; "/Imonospace", [0xFF29]; "/Iniarmenian", [0x053B];
"/Iocyrillic", [0x0401]; "/Iogonek", [0x012E]; "/Iota", [0x0399];
"/Iotaafrican", [0x0196]; "/Iotadieresis", [0x03AA]; "/Iotatonos", [0x038A];
"/Ismall", [0xF769]; "/Istroke", [0x0197]; "/Itilde", [0x0128]; "/Itildebelow",
[0x1E2C]; "/Izhitsacyrillic", [0x0474]; "/Izhitsadblgravecyrillic", [0x0476];
"/J", [0x004A]; "/Jaarmenian", [0x0541]; "/Jcircle", [0x24BF]; "/Jcircumflex",
[0x0134]; "/Jecyrillic", [0x0408]; "/Jheharmenian", [0x054B]; "/Jmonospace",
[0xFF2A]; "/Jsmall", [0xF76A]; "/K", [0x004B]; "/KBsquare", [0x3385];
"/KKsquare", [0x33CD]; "/Kabashkircyrillic", [0x04A0]; "/Kacute", [0x1E30];
"/Kacyrillic", [0x041A]; "/Kadescendercyrillic", [0x049A]; "/Kahookcyrillic",
[0x04C3]; "/Kappa", [0x039A]; "/Kastrokecyrillic", [0x049E];
"/Kaverticalstrokecyrillic", [0x049C]; "/Kcaron", [0x01E8]; "/Kcedilla",
[0x0136]; "/Kcircle", [0x24C0]; "/Kcommaaccent", [0x0136]; "/Kdotbelow",
[0x1E32]; "/Keharmenian", [0x0554]; "/Kenarmenian", [0x053F]; "/Khacyrillic",
[0x0425]; "/Kheicoptic", [0x03E6]; "/Khook", [0x0198]; "/Kjecyrillic", [0x040C];
"/Klinebelow", [0x1E34]; "/Kmonospace", [0xFF2B]; "/Koppacyrillic", [0x0480];
"/Koppagreek", [0x03DE]; "/Ksicyrillic", [0x046E]; "/Ksmall", [0xF76B]; "/L",
[0x004C]; "/LJ", [0x01C7]; "/LL", [0xF6BF]; "/Lacute", [0x0139]; "/Lambda",
[0x039B]; "/Lcaron", [0x013D]; "/Lcedilla", [0x013B]; "/Lcircle", [0x24C1];
"/Lcircumflexbelow", [0x1E3C]; "/Lcommaaccent", [0x013B]; "/Ldot", [0x013F];
"/Ldotaccent", [0x013F]; "/Ldotbelow", [0x1E36]; "/Ldotbelowmacron", [0x1E38];
"/Liwnarmenian", [0x053C]; "/Lj", [0x01C8]; "/Ljecyrillic", [0x0409];
"/Llinebelow", [0x1E3A]; "/Lmonospace", [0xFF2C]; "/Lslash", [0x0141];
"/Lslashsmall", [0xF6F9]; "/Lsmall", [0xF76C]; "/M", [0x004D]; "/MBsquare",
[0x3386]; "/Macron", [0xF6D0]; "/Macronsmall", [0xF7AF]; "/Macute", [0x1E3E];
"/Mcircle", [0x24C2]; "/Mdotaccent", [0x1E40]; "/Mdotbelow", [0x1E42];
"/Menarmenian", [0x0544]; "/Mmonospace", [0xFF2D]; "/Msmall", [0xF76D];
"/Mturned", [0x019C]; "/Mu", [0x039C]; "/N", [0x004E]; "/NJ", [0x01CA];
"/Nacute", [0x0143]; "/Ncaron", [0x0147]; "/Ncedilla", [0x0145]; "/Ncircle",
[0x24C3]; "/Ncircumflexbelow", [0x1E4A]; "/Ncommaaccent", [0x0145];
"/Ndotaccent", [0x1E44]; "/Ndotbelow", [0x1E46]; "/Nhookleft", [0x019D];
"/Nineroman", [0x2168]; "/Nj", [0x01CB]; "/Njecyrillic", [0x040A];
"/Nlinebelow", [0x1E48]; "/Nmonospace", [0xFF2E]; "/Nowarmenian", [0x0546];
"/Nsmall", [0xF76E]; "/Ntilde", [0x00D1]; "/Ntildesmall", [0xF7F1]; "/Nu",
[0x039D]; "/O", [0x004F]; "/OE", [0x0152]; "/OEsmall", [0xF6FA]; "/Oacute",
[0x00D3]; "/Oacutesmall", [0xF7F3]; "/Obarredcyrillic", [0x04E8];
"/Obarreddieresiscyrillic", [0x04EA]; "/Obreve", [0x014E]; "/Ocaron", [0x01D1];
"/Ocenteredtilde", [0x019F]; "/Ocircle", [0x24C4]; "/Ocircumflex", [0x00D4];
"/Ocircumflexacute", [0x1ED0]; "/Ocircumflexdotbelow", [0x1ED8];
"/Ocircumflexgrave", [0x1ED2]; "/Ocircumflexhookabove", [0x1ED4];
"/Ocircumflexsmall", [0xF7F4]; "/Ocircumflextilde", [0x1ED6]; "/Ocyrillic",
[0x041E]; "/Odblacute", [0x0150]; "/Odblgrave", [0x020C]; "/Odieresis",
[0x00D6]; "/Odieresiscyrillic", [0x04E6]; "/Odieresissmall", [0xF7F6];
"/Odotbelow", [0x1ECC]; "/Ogoneksmall", [0xF6FB]; "/Ograve", [0x00D2];
"/Ogravesmall", [0xF7F2]; "/Oharmenian", [0x0555]; "/Ohm", [0x2126];
"/Ohookabove", [0x1ECE]; "/Ohorn", [0x01A0]; "/Ohornacute", [0x1EDA];
"/Ohorndotbelow", [0x1EE2]; "/Ohorngrave", [0x1EDC]; "/Ohornhookabove",
[0x1EDE]; "/Ohorntilde", [0x1EE0]; "/Ohungarumlaut", [0x0150]; "/Oi", [0x01A2];
"/Oinvertedbreve", [0x020E]; "/Omacron", [0x014C]; "/Omacronacute", [0x1E52];
"/Omacrongrave", [0x1E50]; "/Omega", [0x2126]; "/Omegacyrillic", [0x0460];
"/Omegagreek", [0x03A9]; "/Omegaroundcyrillic", [0x047A]; "/Omegatitlocyrillic",
[0x047C]; "/Omegatonos", [0x038F]; "/Omicron", [0x039F]; "/Omicrontonos",
[0x038C]; "/Omonospace", [0xFF2F]; "/Oneroman", [0x2160]; "/Oogonek", [0x01EA];
"/Oogonekmacron", [0x01EC]; "/Oopen", [0x0186]; "/Oslash", [0x00D8];
"/Oslashacute", [0x01FE]; "/Oslashsmall", [0xF7F8]; "/Osmall", [0xF76F];
"/Ostrokeacute", [0x01FE]; "/Otcyrillic", [0x047E]; "/Otilde", [0x00D5];
"/Otildeacute", [0x1E4C]; "/Otildedieresis", [0x1E4E]; "/Otildesmall", [0xF7F5];
"/P", [0x0050]; "/Pacute", [0x1E54]; "/Pcircle", [0x24C5]; "/Pdotaccent",
[0x1E56]; "/Pecyrillic", [0x041F]; "/Peharmenian", [0x054A];
"/Pemiddlehookcyrillic", [0x04A6]; "/Phi", [0x03A6]; "/Phook", [0x01A4]; "/Pi",
[0x03A0]; "/Piwrarmenian", [0x0553]; "/Pmonospace", [0xFF30]; "/Psi", [0x03A8];
"/Psicyrillic", [0x0470]; "/Psmall", [0xF770]; "/Q", [0x0051]; "/Qcircle",
[0x24C6]; "/Qmonospace", [0xFF31]; "/Qsmall", [0xF771]; "/R", [0x0052];
"/Raarmenian", [0x054C]; "/Racute", [0x0154]; "/Rcaron", [0x0158]; "/Rcedilla",
[0x0156]; "/Rcircle", [0x24C7]; "/Rcommaaccent", [0x0156]; "/Rdblgrave",
[0x0210]; "/Rdotaccent", [0x1E58]; "/Rdotbelow", [0x1E5A]; "/Rdotbelowmacron",
[0x1E5C]; "/Reharmenian", [0x0550]; "/Rfraktur", [0x211C]; "/Rho", [0x03A1];
"/Ringsmall", [0xF6FC]; "/Rinvertedbreve", [0x0212]; "/Rlinebelow", [0x1E5E];
"/Rmonospace", [0xFF32]; "/Rsmall", [0xF772]; "/Rsmallinverted", [0x0281];
"/Rsmallinvertedsuperior", [0x02B6]; "/S", [0x0053]; "/SF010000", [0x250C];
"/SF020000", [0x2514]; "/SF030000", [0x2510]; "/SF040000", [0x2518];
"/SF050000", [0x253C]; "/SF060000", [0x252C]; "/SF070000", [0x2534];
"/SF080000", [0x251C]; "/SF090000", [0x2524]; "/SF100000", [0x2500];
"/SF110000", [0x2502]; "/SF190000", [0x2561]; "/SF200000", [0x2562];
"/SF210000", [0x2556]; "/SF220000", [0x2555]; "/SF230000", [0x2563];
"/SF240000", [0x2551]; "/SF250000", [0x2557]; "/SF260000", [0x255D];
"/SF270000", [0x255C]; "/SF280000", [0x255B]; "/SF360000", [0x255E];
"/SF370000", [0x255F]; "/SF380000", [0x255A]; "/SF390000", [0x2554];
"/SF400000", [0x2569]; "/SF410000", [0x2566]; "/SF420000", [0x2560];
"/SF430000", [0x2550]; "/SF440000", [0x256C]; "/SF450000", [0x2567];
"/SF460000", [0x2568]; "/SF470000", [0x2564]; "/SF480000", [0x2565];
"/SF490000", [0x2559]; "/SF500000", [0x2558]; "/SF510000", [0x2552];
"/SF520000", [0x2553]; "/SF530000", [0x256B]; "/SF540000", [0x256A]; "/Sacute",
[0x015A]; "/Sacutedotaccent", [0x1E64]; "/Sampigreek", [0x03E0]; "/Scaron",
[0x0160]; "/Scarondotaccent", [0x1E66]; "/Scaronsmall", [0xF6FD]; "/Scedilla",
[0x015E]; "/Schwa", [0x018F]; "/Schwacyrillic", [0x04D8];
"/Schwadieresiscyrillic", [0x04DA]; "/Scircle", [0x24C8]; "/Scircumflex",
[0x015C]; "/Scommaaccent", [0x0218]; "/Sdotaccent", [0x1E60]; "/Sdotbelow",
[0x1E62]; "/Sdotbelowdotaccent", [0x1E68]; "/Seharmenian", [0x054D];
"/Sevenroman", [0x2166]; "/Shaarmenian", [0x0547]; "/Shacyrillic", [0x0428];
"/Shchacyrillic", [0x0429]; "/Sheicoptic", [0x03E2]; "/Shhacyrillic", [0x04BA];
"/Shimacoptic", [0x03EC]; "/Sigma", [0x03A3]; "/Sixroman", [0x2165];
"/Smonospace", [0xFF33]; "/Softsigncyrillic", [0x042C]; "/Ssmall", [0xF773];
"/Stigmagreek", [0x03DA]; "/T", [0x0054]; "/Tau", [0x03A4]; "/Tbar", [0x0166];
"/Tcaron", [0x0164]; "/Tcedilla", [0x0162]; "/Tcircle", [0x24C9];
"/Tcircumflexbelow", [0x1E70]; "/Tcommaaccent", [0x0162]; "/Tdotaccent",
[0x1E6A]; "/Tdotbelow", [0x1E6C]; "/Tecyrillic", [0x0422];
"/Tedescendercyrillic", [0x04AC]; "/Tenroman", [0x2169]; "/Tetsecyrillic",
[0x04B4]; "/Theta", [0x0398]; "/Thook", [0x01AC]; "/Thorn", [0x00DE];
"/Thornsmall", [0xF7FE]; "/Threeroman", [0x2162]; "/Tildesmall", [0xF6FE];
"/Tiwnarmenian", [0x054F]; "/Tlinebelow", [0x1E6E]; "/Tmonospace", [0xFF34];
"/Toarmenian", [0x0539]; "/Tonefive", [0x01BC]; "/Tonesix", [0x0184];
"/Tonetwo", [0x01A7]; "/Tretroflexhook", [0x01AE]; "/Tsecyrillic", [0x0426];
"/Tshecyrillic", [0x040B]; "/Tsmall", [0xF774]; "/Twelveroman", [0x216B];
"/Tworoman", [0x2161]; "/U", [0x0055]; "/Uacute", [0x00DA]; "/Uacutesmall",
[0xF7FA]; "/Ubreve", [0x016C]; "/Ucaron", [0x01D3]; "/Ucircle", [0x24CA];
"/Ucircumflex", [0x00DB]; "/Ucircumflexbelow", [0x1E76]; "/Ucircumflexsmall",
[0xF7FB]; "/Ucyrillic", [0x0423]; "/Udblacute", [0x0170]; "/Udblgrave",
[0x0214]; "/Udieresis", [0x00DC]; "/Udieresisacute", [0x01D7];
"/Udieresisbelow", [0x1E72]; "/Udieresiscaron", [0x01D9]; "/Udieresiscyrillic",
[0x04F0]; "/Udieresisgrave", [0x01DB]; "/Udieresismacron", [0x01D5];
"/Udieresissmall", [0xF7FC]; "/Udotbelow", [0x1EE4]; "/Ugrave", [0x00D9];
"/Ugravesmall", [0xF7F9]; "/Uhookabove", [0x1EE6]; "/Uhorn", [0x01AF];
"/Uhornacute", [0x1EE8]; "/Uhorndotbelow", [0x1EF0]; "/Uhorngrave", [0x1EEA];
"/Uhornhookabove", [0x1EEC]; "/Uhorntilde", [0x1EEE]; "/Uhungarumlaut",
[0x0170]; "/Uhungarumlautcyrillic", [0x04F2]; "/Uinvertedbreve", [0x0216];
"/Ukcyrillic", [0x0478]; "/Umacron", [0x016A]; "/Umacroncyrillic", [0x04EE];
"/Umacrondieresis", [0x1E7A]; "/Umonospace", [0xFF35]; "/Uogonek", [0x0172];
"/Upsilon", [0x03A5]; "/Upsilon1", [0x03D2]; "/Upsilonacutehooksymbolgreek",
[0x03D3]; "/Upsilonafrican", [0x01B1]; "/Upsilondieresis", [0x03AB];
"/Upsilondieresishooksymbolgreek", [0x03D4]; "/Upsilonhooksymbol", [0x03D2];
"/Upsilontonos", [0x038E]; "/Uring", [0x016E]; "/Ushortcyrillic", [0x040E];
"/Usmall", [0xF775]; "/Ustraightcyrillic", [0x04AE]; "/Ustraightstrokecyrillic",
[0x04B0]; "/Utilde", [0x0168]; "/Utildeacute", [0x1E78]; "/Utildebelow",
[0x1E74]; "/V", [0x0056]; "/Vcircle", [0x24CB]; "/Vdotbelow", [0x1E7E];
"/Vecyrillic", [0x0412]; "/Vewarmenian", [0x054E]; "/Vhook", [0x01B2];
"/Vmonospace", [0xFF36]; "/Voarmenian", [0x0548]; "/Vsmall", [0xF776];
"/Vtilde", [0x1E7C]; "/W", [0x0057]; "/Wacute", [0x1E82]; "/Wcircle", [0x24CC];
"/Wcircumflex", [0x0174]; "/Wdieresis", [0x1E84]; "/Wdotaccent", [0x1E86];
"/Wdotbelow", [0x1E88]; "/Wgrave", [0x1E80]; "/Wmonospace", [0xFF37]; "/Wsmall",
[0xF777]; "/X", [0x0058]; "/Xcircle", [0x24CD]; "/Xdieresis", [0x1E8C];
"/Xdotaccent", [0x1E8A]; "/Xeharmenian", [0x053D]; "/Xi", [0x039E];
"/Xmonospace", [0xFF38]; "/Xsmall", [0xF778]; "/Y", [0x0059]; "/Yacute",
[0x00DD]; "/Yacutesmall", [0xF7FD]; "/Yatcyrillic", [0x0462]; "/Ycircle",
[0x24CE]; "/Ycircumflex", [0x0176]; "/Ydieresis", [0x0178]; "/Ydieresissmall",
[0xF7FF]; "/Ydotaccent", [0x1E8E]; "/Ydotbelow", [0x1EF4]; "/Yericyrillic",
[0x042B]; "/Yerudieresiscyrillic", [0x04F8]; "/Ygrave", [0x1EF2]; "/Yhook",
[0x01B3]; "/Yhookabove", [0x1EF6]; "/Yiarmenian", [0x0545]; "/Yicyrillic",
[0x0407]; "/Yiwnarmenian", [0x0552]; "/Ymonospace", [0xFF39]; "/Ysmall",
[0xF779]; "/Ytilde", [0x1EF8]; "/Yusbigcyrillic", [0x046A];
"/Yusbigiotifiedcyrillic", [0x046C]; "/Yuslittlecyrillic", [0x0466];
"/Yuslittleiotifiedcyrillic", [0x0468]; "/Z", [0x005A]; "/Zaarmenian", [0x0536];
"/Zacute", [0x0179]; "/Zcaron", [0x017D]; "/Zcaronsmall", [0xF6FF]; "/Zcircle",
[0x24CF]; "/Zcircumflex", [0x1E90]; "/Zdot", [0x017B]; "/Zdotaccent", [0x017B];
"/Zdotbelow", [0x1E92]; "/Zecyrillic", [0x0417]; "/Zedescendercyrillic",
[0x0498]; "/Zedieresiscyrillic", [0x04DE]; "/Zeta", [0x0396]; "/Zhearmenian",
[0x053A]; "/Zhebrevecyrillic", [0x04C1]; "/Zhecyrillic", [0x0416];
"/Zhedescendercyrillic", [0x0496]; "/Zhedieresiscyrillic", [0x04DC];
"/Zlinebelow", [0x1E94]; "/Zmonospace", [0xFF3A]; "/Zsmall", [0xF77A];
"/Zstroke", [0x01B5]; "/a", [0x0061]; "/aabengali", [0x0986]; "/aacute",
[0x00E1]; "/aadeva", [0x0906]; "/aagujarati", [0x0A86]; "/aagurmukhi", [0x0A06];
"/aamatragurmukhi", [0x0A3E]; "/aarusquare", [0x3303]; "/aavowelsignbengali",
[0x09BE]; "/aavowelsigndeva", [0x093E]; "/aavowelsigngujarati", [0x0ABE];
"/abbreviationmarkarmenian", [0x055F]; "/abbreviationsigndeva", [0x0970];
"/abengali", [0x0985]; "/abopomofo", [0x311A]; "/abreve", [0x0103];
"/abreveacute", [0x1EAF]; "/abrevecyrillic", [0x04D1]; "/abrevedotbelow",
[0x1EB7]; "/abrevegrave", [0x1EB1]; "/abrevehookabove", [0x1EB3];
"/abrevetilde", [0x1EB5]; "/acaron", [0x01CE]; "/acircle", [0x24D0];
"/acircumflex", [0x00E2]; "/acircumflexacute", [0x1EA5]; "/acircumflexdotbelow",
[0x1EAD]; "/acircumflexgrave", [0x1EA7]; "/acircumflexhookabove", [0x1EA9];
"/acircumflextilde", [0x1EAB]; "/acute", [0x00B4]; "/acutebelowcmb", [0x0317];
"/acutecmb", [0x0301]; "/acutecomb", [0x0301]; "/acutedeva", [0x0954];
"/acutelowmod", [0x02CF]; "/acutetonecmb", [0x0341]; "/acyrillic", [0x0430];
"/adblgrave", [0x0201]; "/addakgurmukhi", [0x0A71]; "/adeva", [0x0905];
"/adieresis", [0x00E4]; "/adieresiscyrillic", [0x04D3]; "/adieresismacron",
[0x01DF]; "/adotbelow", [0x1EA1]; "/adotmacron", [0x01E1]; "/ae", [0x00E6];
"/aeacute", [0x01FD]; "/aekorean", [0x3150]; "/aemacron", [0x01E3];
"/afii00208", [0x2015]; "/afii08941", [0x20A4]; "/afii10017", [0x0410];
"/afii10018", [0x0411]; "/afii10019", [0x0412]; "/afii10020", [0x0413];
"/afii10021", [0x0414]; "/afii10022", [0x0415]; "/afii10023", [0x0401];
"/afii10024", [0x0416]; "/afii10025", [0x0417]; "/afii10026", [0x0418];
"/afii10027", [0x0419]; "/afii10028", [0x041A]; "/afii10029", [0x041B];
"/afii10030", [0x041C]; "/afii10031", [0x041D]; "/afii10032", [0x041E];
"/afii10033", [0x041F]; "/afii10034", [0x0420]; "/afii10035", [0x0421];
"/afii10036", [0x0422]; "/afii10037", [0x0423]; "/afii10038", [0x0424];
"/afii10039", [0x0425]; "/afii10040", [0x0426]; "/afii10041", [0x0427];
"/afii10042", [0x0428]; "/afii10043", [0x0429]; "/afii10044", [0x042A];
"/afii10045", [0x042B]; "/afii10046", [0x042C]; "/afii10047", [0x042D];
"/afii10048", [0x042E]; "/afii10049", [0x042F]; "/afii10050", [0x0490];
"/afii10051", [0x0402]; "/afii10052", [0x0403]; "/afii10053", [0x0404];
"/afii10054", [0x0405]; "/afii10055", [0x0406]; "/afii10056", [0x0407];
"/afii10057", [0x0408]; "/afii10058", [0x0409]; "/afii10059", [0x040A];
"/afii10060", [0x040B]; "/afii10061", [0x040C]; "/afii10062", [0x040E];
"/afii10063", [0xF6C4]; "/afii10064", [0xF6C5]; "/afii10065", [0x0430];
"/afii10066", [0x0431]; "/afii10067", [0x0432]; "/afii10068", [0x0433];
"/afii10069", [0x0434]; "/afii10070", [0x0435]; "/afii10071", [0x0451];
"/afii10072", [0x0436]; "/afii10073", [0x0437]; "/afii10074", [0x0438];
"/afii10075", [0x0439]; "/afii10076", [0x043A]; "/afii10077", [0x043B];
"/afii10078", [0x043C]; "/afii10079", [0x043D]; "/afii10080", [0x043E];
"/afii10081", [0x043F]; "/afii10082", [0x0440]; "/afii10083", [0x0441];
"/afii10084", [0x0442]; "/afii10085", [0x0443]; "/afii10086", [0x0444];
"/afii10087", [0x0445]; "/afii10088", [0x0446]; "/afii10089", [0x0447];
"/afii10090", [0x0448]; "/afii10091", [0x0449]; "/afii10092", [0x044A];
"/afii10093", [0x044B]; "/afii10094", [0x044C]; "/afii10095", [0x044D];
"/afii10096", [0x044E]; "/afii10097", [0x044F]; "/afii10098", [0x0491];
"/afii10099", [0x0452]; "/afii10100", [0x0453]; "/afii10101", [0x0454];
"/afii10102", [0x0455]; "/afii10103", [0x0456]; "/afii10104", [0x0457];
"/afii10105", [0x0458]; "/afii10106", [0x0459]; "/afii10107", [0x045A];
"/afii10108", [0x045B]; "/afii10109", [0x045C]; "/afii10110", [0x045E];
"/afii10145", [0x040F]; "/afii10146", [0x0462]; "/afii10147", [0x0472];
"/afii10148", [0x0474]; "/afii10192", [0xF6C6]; "/afii10193", [0x045F];
"/afii10194", [0x0463]; "/afii10195", [0x0473]; "/afii10196", [0x0475];
"/afii10831", [0xF6C7]; "/afii10832", [0xF6C8]; "/afii10846", [0x04D9];
"/afii299", [0x200E]; "/afii300", [0x200F]; "/afii301", [0x200D]; "/afii57381",
[0x066A]; "/afii57388", [0x060C]; "/afii57392", [0x0660]; "/afii57393",
[0x0661]; "/afii57394", [0x0662]; "/afii57395", [0x0663]; "/afii57396",
[0x0664]; "/afii57397", [0x0665]; "/afii57398", [0x0666]; "/afii57399",
[0x0667]; "/afii57400", [0x0668]; "/afii57401", [0x0669]; "/afii57403",
[0x061B]; "/afii57407", [0x061F]; "/afii57409", [0x0621]; "/afii57410",
[0x0622]; "/afii57411", [0x0623]; "/afii57412", [0x0624]; "/afii57413",
[0x0625]; "/afii57414", [0x0626]; "/afii57415", [0x0627]; "/afii57416",
[0x0628]; "/afii57417", [0x0629]; "/afii57418", [0x062A]; "/afii57419",
[0x062B]; "/afii57420", [0x062C]; "/afii57421", [0x062D]; "/afii57422",
[0x062E]; "/afii57423", [0x062F]; "/afii57424", [0x0630]; "/afii57425",
[0x0631]; "/afii57426", [0x0632]; "/afii57427", [0x0633]; "/afii57428",
[0x0634]; "/afii57429", [0x0635]; "/afii57430", [0x0636]; "/afii57431",
[0x0637]; "/afii57432", [0x0638]; "/afii57433", [0x0639]; "/afii57434",
[0x063A]; "/afii57440", [0x0640]; "/afii57441", [0x0641]; "/afii57442",
[0x0642]; "/afii57443", [0x0643]; "/afii57444", [0x0644]; "/afii57445",
[0x0645]; "/afii57446", [0x0646]; "/afii57448", [0x0648]; "/afii57449",
[0x0649]; "/afii57450", [0x064A]; "/afii57451", [0x064B]; "/afii57452",
[0x064C]; "/afii57453", [0x064D]; "/afii57454", [0x064E]; "/afii57455",
[0x064F]; "/afii57456", [0x0650]; "/afii57457", [0x0651]; "/afii57458",
[0x0652]; "/afii57470", [0x0647]; "/afii57505", [0x06A4]; "/afii57506",
[0x067E]; "/afii57507", [0x0686]; "/afii57508", [0x0698]; "/afii57509",
[0x06AF]; "/afii57511", [0x0679]; "/afii57512", [0x0688]; "/afii57513",
[0x0691]; "/afii57514", [0x06BA]; "/afii57519", [0x06D2]; "/afii57534",
[0x06D5]; "/afii57636", [0x20AA]; "/afii57645", [0x05BE]; "/afii57658",
[0x05C3]; "/afii57664", [0x05D0]; "/afii57665", [0x05D1]; "/afii57666",
[0x05D2]; "/afii57667", [0x05D3]; "/afii57668", [0x05D4]; "/afii57669",
[0x05D5]; "/afii57670", [0x05D6]; "/afii57671", [0x05D7]; "/afii57672",
[0x05D8]; "/afii57673", [0x05D9]; "/afii57674", [0x05DA]; "/afii57675",
[0x05DB]; "/afii57676", [0x05DC]; "/afii57677", [0x05DD]; "/afii57678",
[0x05DE]; "/afii57679", [0x05DF]; "/afii57680", [0x05E0]; "/afii57681",
[0x05E1]; "/afii57682", [0x05E2]; "/afii57683", [0x05E3]; "/afii57684",
[0x05E4]; "/afii57685", [0x05E5]; "/afii57686", [0x05E6]; "/afii57687",
[0x05E7]; "/afii57688", [0x05E8]; "/afii57689", [0x05E9]; "/afii57690",
[0x05EA]; "/afii57694", [0xFB2A]; "/afii57695", [0xFB2B]; "/afii57700",
[0xFB4B]; "/afii57705", [0xFB1F]; "/afii57716", [0x05F0]; "/afii57717",
[0x05F1]; "/afii57718", [0x05F2]; "/afii57723", [0xFB35]; "/afii57793",
[0x05B4]; "/afii57794", [0x05B5]; "/afii57795", [0x05B6]; "/afii57796",
[0x05BB]; "/afii57797", [0x05B8]; "/afii57798", [0x05B7]; "/afii57799",
[0x05B0]; "/afii57800", [0x05B2]; "/afii57801", [0x05B1]; "/afii57802",
[0x05B3]; "/afii57803", [0x05C2]; "/afii57804", [0x05C1]; "/afii57806",
[0x05B9]; "/afii57807", [0x05BC]; "/afii57839", [0x05BD]; "/afii57841",
[0x05BF]; "/afii57842", [0x05C0]; "/afii57929", [0x02BC]; "/afii61248",
[0x2105]; "/afii61289", [0x2113]; "/afii61352", [0x2116]; "/afii61573",
[0x202C]; "/afii61574", [0x202D]; "/afii61575", [0x202E]; "/afii61664",
[0x200C]; "/afii63167", [0x066D]; "/afii64937", [0x02BD]; "/agrave", [0x00E0];
"/agujarati", [0x0A85]; "/agurmukhi", [0x0A05]; "/ahiragana", [0x3042];
"/ahookabove", [0x1EA3]; "/aibengali", [0x0990]; "/aibopomofo", [0x311E];
"/aideva", [0x0910]; "/aiecyrillic", [0x04D5]; "/aigujarati", [0x0A90];
"/aigurmukhi", [0x0A10]; "/aimatragurmukhi", [0x0A48]; "/ainarabic", [0x0639];
"/ainfinalarabic", [0xFECA]; "/aininitialarabic", [0xFECB]; "/ainmedialarabic",
[0xFECC]; "/ainvertedbreve", [0x0203]; "/aivowelsignbengali", [0x09C8];
"/aivowelsigndeva", [0x0948]; "/aivowelsigngujarati", [0x0AC8]; "/akatakana",
[0x30A2]; "/akatakanahalfwidth", [0xFF71]; "/akorean", [0x314F]; "/alef",
[0x05D0]; "/alefarabic", [0x0627]; "/alefdageshhebrew", [0xFB30];
"/aleffinalarabic", [0xFE8E]; "/alefhamzaabovearabic", [0x0623];
"/alefhamzaabovefinalarabic", [0xFE84]; "/alefhamzabelowarabic", [0x0625];
"/alefhamzabelowfinalarabic", [0xFE88]; "/alefhebrew", [0x05D0];
"/aleflamedhebrew", [0xFB4F]; "/alefmaddaabovearabic", [0x0622];
"/alefmaddaabovefinalarabic", [0xFE82]; "/alefmaksuraarabic", [0x0649];
"/alefmaksurafinalarabic", [0xFEF0]; "/alefmaksurainitialarabic", [0xFEF3];
"/alefmaksuramedialarabic", [0xFEF4]; "/alefpatahhebrew", [0xFB2E];
"/alefqamatshebrew", [0xFB2F]; "/aleph", [0x2135]; "/allequal", [0x224C];
"/alpha", [0x03B1]; "/alphatonos", [0x03AC]; "/amacron", [0x0101];
"/amonospace", [0xFF41]; "/ampersand", [0x0026]; "/ampersandmonospace",
[0xFF06]; "/ampersandsmall", [0xF726]; "/amsquare", [0x33C2]; "/anbopomofo",
[0x3122]; "/angbopomofo", [0x3124]; "/angkhankhuthai", [0x0E5A]; "/angle",
[0x2220]; "/anglebracketleft", [0x3008]; "/anglebracketleftvertical", [0xFE3F];
"/anglebracketright", [0x3009]; "/anglebracketrightvertical", [0xFE40];
"/angleleft", [0x2329]; "/angleright", [0x232A]; "/angstrom", [0x212B];
"/anoteleia", [0x0387]; "/anudattadeva", [0x0952]; "/anusvarabengali", [0x0982];
"/anusvaradeva", [0x0902]; "/anusvaragujarati", [0x0A82]; "/aogonek", [0x0105];
"/apaatosquare", [0x3300]; "/aparen", [0x249C]; "/apostrophearmenian", [0x055A];
"/apostrophemod", [0x02BC]; "/apple", [0xF8FF]; "/approaches", [0x2250];
"/approxequal", [0x2248]; "/approxequalorimage", [0x2252];
"/approximatelyequal", [0x2245]; "/araeaekorean", [0x318E]; "/araeakorean",
[0x318D]; "/arc", [0x2312]; "/arighthalfring", [0x1E9A]; "/aring", [0x00E5];
"/aringacute", [0x01FB]; "/aringbelow", [0x1E01]; "/arrowboth", [0x2194];
"/arrowdashdown", [0x21E3]; "/arrowdashleft", [0x21E0]; "/arrowdashright",
[0x21E2]; "/arrowdashup", [0x21E1]; "/arrowdblboth", [0x21D4]; "/arrowdbldown",
[0x21D3]; "/arrowdblleft", [0x21D0]; "/arrowdblright", [0x21D2]; "/arrowdblup",
[0x21D1]; "/arrowdown", [0x2193]; "/arrowdownleft", [0x2199]; "/arrowdownright",
[0x2198]; "/arrowdownwhite", [0x21E9]; "/arrowheaddownmod", [0x02C5];
"/arrowheadleftmod", [0x02C2]; "/arrowheadrightmod", [0x02C3];
"/arrowheadupmod", [0x02C4]; "/arrowhorizex", [0xF8E7]; "/arrowleft", [0x2190];
"/arrowleftdbl", [0x21D0]; "/arrowleftdblstroke", [0x21CD];
"/arrowleftoverright", [0x21C6]; "/arrowleftwhite", [0x21E6]; "/arrowright",
[0x2192]; "/arrowrightdblstroke", [0x21CF]; "/arrowrightheavy", [0x279E];
"/arrowrightoverleft", [0x21C4]; "/arrowrightwhite", [0x21E8]; "/arrowtableft",
[0x21E4]; "/arrowtabright", [0x21E5]; "/arrowup", [0x2191]; "/arrowupdn",
[0x2195]; "/arrowupdnbse", [0x21A8]; "/arrowupdownbase", [0x21A8];
"/arrowupleft", [0x2196]; "/arrowupleftofdown", [0x21C5]; "/arrowupright",
[0x2197]; "/arrowupwhite", [0x21E7]; "/arrowvertex", [0xF8E6]; "/asciicircum",
[0x005E]; "/asciicircummonospace", [0xFF3E]; "/asciitilde", [0x007E];
"/asciitildemonospace", [0xFF5E]; "/ascript", [0x0251]; "/ascriptturned",
[0x0252]; "/asmallhiragana", [0x3041]; "/asmallkatakana", [0x30A1];
"/asmallkatakanahalfwidth", [0xFF67]; "/asterisk", [0x002A];
"/asteriskaltonearabic", [0x066D]; "/asteriskarabic", [0x066D]; "/asteriskmath",
[0x2217]; "/asteriskmonospace", [0xFF0A]; "/asterisksmall", [0xFE61];
"/asterism", [0x2042]; "/asuperior", [0xF6E9]; "/asymptoticallyequal", [0x2243];
"/at", [0x0040]; "/atilde", [0x00E3]; "/atmonospace", [0xFF20]; "/atsmall",
[0xFE6B]; "/aturned", [0x0250]; "/aubengali", [0x0994]; "/aubopomofo", [0x3120];
"/audeva", [0x0914]; "/augujarati", [0x0A94]; "/augurmukhi", [0x0A14];
"/aulengthmarkbengali", [0x09D7]; "/aumatragurmukhi", [0x0A4C];
"/auvowelsignbengali", [0x09CC]; "/auvowelsigndeva", [0x094C];
"/auvowelsigngujarati", [0x0ACC]; "/avagrahadeva", [0x093D]; "/aybarmenian",
[0x0561]; "/ayin", [0x05E2]; "/ayinaltonehebrew", [0xFB20]; "/ayinhebrew",
[0x05E2]; "/b", [0x0062]; "/babengali", [0x09AC]; "/backslash", [0x005C];
"/backslashmonospace", [0xFF3C]; "/badeva", [0x092C]; "/bagujarati", [0x0AAC];
"/bagurmukhi", [0x0A2C]; "/bahiragana", [0x3070]; "/bahtthai", [0x0E3F];
"/bakatakana", [0x30D0]; "/bar", [0x007C]; "/barmonospace", [0xFF5C];
"/bbopomofo", [0x3105]; "/bcircle", [0x24D1]; "/bdotaccent", [0x1E03];
"/bdotbelow", [0x1E05]; "/beamedsixteenthnotes", [0x266C]; "/because", [0x2235];
"/becyrillic", [0x0431]; "/beharabic", [0x0628]; "/behfinalarabic", [0xFE90];
"/behinitialarabic", [0xFE91]; "/behiragana", [0x3079]; "/behmedialarabic",
[0xFE92]; "/behmeeminitialarabic", [0xFC9F]; "/behmeemisolatedarabic", [0xFC08];
"/behnoonfinalarabic", [0xFC6D]; "/bekatakana", [0x30D9]; "/benarmenian",
[0x0562]; "/bet", [0x05D1]; "/beta", [0x03B2]; "/betasymbolgreek", [0x03D0];
"/betdagesh", [0xFB31]; "/betdageshhebrew", [0xFB31]; "/bethebrew", [0x05D1];
"/betrafehebrew", [0xFB4C]; "/bhabengali", [0x09AD]; "/bhadeva", [0x092D];
"/bhagujarati", [0x0AAD]; "/bhagurmukhi", [0x0A2D]; "/bhook", [0x0253];
"/bihiragana", [0x3073]; "/bikatakana", [0x30D3]; "/bilabialclick", [0x0298];
"/bindigurmukhi", [0x0A02]; "/birusquare", [0x3331]; "/blackcircle", [0x25CF];
"/blackdiamond", [0x25C6]; "/blackdownpointingtriangle", [0x25BC];
"/blackleftpointingpointer", [0x25C4]; "/blackleftpointingtriangle", [0x25C0];
"/blacklenticularbracketleft", [0x3010]; "/blacklenticularbracketleftvertical",
[0xFE3B]; "/blacklenticularbracketright", [0x3011];
"/blacklenticularbracketrightvertical", [0xFE3C]; "/blacklowerlefttriangle",
[0x25E3]; "/blacklowerrighttriangle", [0x25E2]; "/blackrectangle", [0x25AC];
"/blackrightpointingpointer", [0x25BA]; "/blackrightpointingtriangle", [0x25B6];
"/blacksmallsquare", [0x25AA]; "/blacksmilingface", [0x263B]; "/blacksquare",
[0x25A0]; "/blackstar", [0x2605]; "/blackupperlefttriangle", [0x25E4];
"/blackupperrighttriangle", [0x25E5]; "/blackuppointingsmalltriangle", [0x25B4];
"/blackuppointingtriangle", [0x25B2]; "/blank", [0x2423]; "/blinebelow",
[0x1E07]; "/block", [0x2588]; "/bmonospace", [0xFF42]; "/bobaimaithai",
[0x0E1A]; "/bohiragana", [0x307C]; "/bokatakana", [0x30DC]; "/bparen", [0x249D];
"/bqsquare", [0x33C3]; "/braceex", [0xF8F4]; "/braceleft", [0x007B];
"/braceleftbt", [0xF8F3]; "/braceleftmid", [0xF8F2]; "/braceleftmonospace",
[0xFF5B]; "/braceleftsmall", [0xFE5B]; "/bracelefttp", [0xF8F1];
"/braceleftvertical", [0xFE37]; "/braceright", [0x007D]; "/bracerightbt",
[0xF8FE]; "/bracerightmid", [0xF8FD]; "/bracerightmonospace", [0xFF5D];
"/bracerightsmall", [0xFE5C]; "/bracerighttp", [0xF8FC]; "/bracerightvertical",
[0xFE38]; "/bracketleft", [0x005B]; "/bracketleftbt", [0xF8F0];
"/bracketleftex", [0xF8EF]; "/bracketleftmonospace", [0xFF3B]; "/bracketlefttp",
[0xF8EE]; "/bracketright", [0x005D]; "/bracketrightbt", [0xF8FB];
"/bracketrightex", [0xF8FA]; "/bracketrightmonospace", [0xFF3D];
"/bracketrighttp", [0xF8F9]; "/breve", [0x02D8]; "/brevebelowcmb", [0x032E];
"/brevecmb", [0x0306]; "/breveinvertedbelowcmb", [0x032F]; "/breveinvertedcmb",
[0x0311]; "/breveinverteddoublecmb", [0x0361]; "/bridgebelowcmb", [0x032A];
"/bridgeinvertedbelowcmb", [0x033A]; "/brokenbar", [0x00A6]; "/bstroke",
[0x0180]; "/bsuperior", [0xF6EA]; "/btopbar", [0x0183]; "/buhiragana", [0x3076];
"/bukatakana", [0x30D6]; "/bullet", [0x2022]; "/bulletinverse", [0x25D8];
"/bulletoperator", [0x2219]; "/bullseye", [0x25CE]; "/c", [0x0063];
"/caarmenian", [0x056E]; "/cabengali", [0x099A]; "/cacute", [0x0107]; "/cadeva",
[0x091A]; "/cagujarati", [0x0A9A]; "/cagurmukhi", [0x0A1A]; "/calsquare",
[0x3388]; "/candrabindubengali", [0x0981]; "/candrabinducmb", [0x0310];
"/candrabindudeva", [0x0901]; "/candrabindugujarati", [0x0A81]; "/capslock",
[0x21EA]; "/careof", [0x2105]; "/caron", [0x02C7]; "/caronbelowcmb", [0x032C];
"/caroncmb", [0x030C]; "/carriagereturn", [0x21B5]; "/cbopomofo", [0x3118];
"/ccaron", [0x010D]; "/ccedilla", [0x00E7]; "/ccedillaacute", [0x1E09];
"/ccircle", [0x24D2]; "/ccircumflex", [0x0109]; "/ccurl", [0x0255]; "/cdot",
[0x010B]; "/cdotaccent", [0x010B]; "/cdsquare", [0x33C5]; "/cedilla", [0x00B8];
"/cedillacmb", [0x0327]; "/cent", [0x00A2]; "/centigrade", [0x2103];
"/centinferior", [0xF6DF]; "/centmonospace", [0xFFE0]; "/centoldstyle",
[0xF7A2]; "/centsuperior", [0xF6E0]; "/chaarmenian", [0x0579]; "/chabengali",
[0x099B]; "/chadeva", [0x091B]; "/chagujarati", [0x0A9B]; "/chagurmukhi",
[0x0A1B]; "/chbopomofo", [0x3114]; "/cheabkhasiancyrillic", [0x04BD];
"/checkmark", [0x2713]; "/checyrillic", [0x0447];
"/chedescenderabkhasiancyrillic", [0x04BF]; "/chedescendercyrillic", [0x04B7];
"/chedieresiscyrillic", [0x04F5]; "/cheharmenian", [0x0573];
"/chekhakassiancyrillic", [0x04CC]; "/cheverticalstrokecyrillic", [0x04B9];
"/chi", [0x03C7]; "/chieuchacirclekorean", [0x3277]; "/chieuchaparenkorean",
[0x3217]; "/chieuchcirclekorean", [0x3269]; "/chieuchkorean", [0x314A];
"/chieuchparenkorean", [0x3209]; "/chochangthai", [0x0E0A]; "/chochanthai",
[0x0E08]; "/chochingthai", [0x0E09]; "/chochoethai", [0x0E0C]; "/chook",
[0x0188]; "/cieucacirclekorean", [0x3276]; "/cieucaparenkorean", [0x3216];
"/cieuccirclekorean", [0x3268]; "/cieuckorean", [0x3148]; "/cieucparenkorean",
[0x3208]; "/cieucuparenkorean", [0x321C]; "/circle", [0x25CB];
"/circlemultiply", [0x2297]; "/circleot", [0x2299]; "/circleplus", [0x2295];
"/circlepostalmark", [0x3036]; "/circlewithlefthalfblack", [0x25D0];
"/circlewithrighthalfblack", [0x25D1]; "/circumflex", [0x02C6];
"/circumflexbelowcmb", [0x032D]; "/circumflexcmb", [0x0302]; "/clear", [0x2327];
"/clickalveolar", [0x01C2]; "/clickdental", [0x01C0]; "/clicklateral", [0x01C1];
"/clickretroflex", [0x01C3]; "/club", [0x2663]; "/clubsuitblack", [0x2663];
"/clubsuitwhite", [0x2667]; "/cmcubedsquare", [0x33A4]; "/cmonospace", [0xFF43];
"/cmsquaredsquare", [0x33A0]; "/coarmenian", [0x0581]; "/colon", [0x003A];
"/colonmonetary", [0x20A1]; "/colonmonospace", [0xFF1A]; "/colonsign", [0x20A1];
"/colonsmall", [0xFE55]; "/colontriangularhalfmod", [0x02D1];
"/colontriangularmod", [0x02D0]; "/comma", [0x002C]; "/commaabovecmb", [0x0313];
"/commaaboverightcmb", [0x0315]; "/commaaccent", [0xF6C3]; "/commaarabic",
[0x060C]; "/commaarmenian", [0x055D]; "/commainferior", [0xF6E1];
"/commamonospace", [0xFF0C]; "/commareversedabovecmb", [0x0314];
"/commareversedmod", [0x02BD]; "/commasmall", [0xFE50]; "/commasuperior",
[0xF6E2]; "/commaturnedabovecmb", [0x0312]; "/commaturnedmod", [0x02BB];
"/compass", [0x263C]; "/congruent", [0x2245]; "/contourintegral", [0x222E];
"/control", [0x2303]; "/controlACK", [0x0006]; "/controlBEL", [0x0007];
"/controlBS", [0x0008]; "/controlCAN", [0x0018]; "/controlCR", [0x000D];
"/controlDC1", [0x0011]; "/controlDC2", [0x0012]; "/controlDC3", [0x0013];
"/controlDC4", [0x0014]; "/controlDEL", [0x007F]; "/controlDLE", [0x0010];
"/controlEM", [0x0019]; "/controlENQ", [0x0005]; "/controlEOT", [0x0004];
"/controlESC", [0x001B]; "/controlETB", [0x0017]; "/controlETX", [0x0003];
"/controlFF", [0x000C]; "/controlFS", [0x001C]; "/controlGS", [0x001D];
"/controlHT", [0x0009]; "/controlLF", [0x000A]; "/controlNAK", [0x0015];
"/controlRS", [0x001E]; "/controlSI", [0x000F]; "/controlSO", [0x000E];
"/controlSOT", [0x0002]; "/controlSTX", [0x0001]; "/controlSUB", [0x001A];
"/controlSYN", [0x0016]; "/controlUS", [0x001F]; "/controlVT", [0x000B];
"/copyright", [0x00A9]; "/copyrightsans", [0xF8E9]; "/copyrightserif", [0xF6D9];
"/cornerbracketleft", [0x300C]; "/cornerbracketlefthalfwidth", [0xFF62];
"/cornerbracketleftvertical", [0xFE41]; "/cornerbracketright", [0x300D];
"/cornerbracketrighthalfwidth", [0xFF63]; "/cornerbracketrightvertical",
[0xFE42]; "/corporationsquare", [0x337F]; "/cosquare", [0x33C7];
"/coverkgsquare", [0x33C6]; "/cparen", [0x249E]; "/cruzeiro", [0x20A2];
"/cstretched", [0x0297]; "/curlyand", [0x22CF]; "/curlyor", [0x22CE];
"/currency", [0x00A4]; "/cyrBreve", [0xF6D1]; "/cyrFlex", [0xF6D2]; "/cyrbreve",
[0xF6D4]; "/cyrflex", [0xF6D5]; "/d", [0x0064]; "/daarmenian", [0x0564];
"/dabengali", [0x09A6]; "/dadarabic", [0x0636]; "/dadeva", [0x0926];
"/dadfinalarabic", [0xFEBE]; "/dadinitialarabic", [0xFEBF]; "/dadmedialarabic",
[0xFEC0]; "/dagesh", [0x05BC]; "/dageshhebrew", [0x05BC]; "/dagger", [0x2020];
"/daggerdbl", [0x2021]; "/dagujarati", [0x0AA6]; "/dagurmukhi", [0x0A26];
"/dahiragana", [0x3060]; "/dakatakana", [0x30C0]; "/dalarabic", [0x062F];
"/dalet", [0x05D3]; "/daletdagesh", [0xFB33]; "/daletdageshhebrew", [0xFB33];
"/dalethatafpatah", [0x05D3; 0x05B2]; "/dalethatafpatahhebrew", [0x05D3;
0x05B2]; "/dalethatafsegol", [0x05D3; 0x05B1]; "/dalethatafsegolhebrew",
[0x05D3; 0x05B1]; "/dalethebrew", [0x05D3]; "/dalethiriq", [0x05D3; 0x05B4];
"/dalethiriqhebrew", [0x05D3; 0x05B4]; "/daletholam", [0x05D3; 0x05B9];
"/daletholamhebrew", [0x05D3; 0x05B9]; "/daletpatah", [0x05D3; 0x05B7];
"/daletpatahhebrew", [0x05D3; 0x05B7]; "/daletqamats", [0x05D3; 0x05B8];
"/daletqamatshebrew", [0x05D3; 0x05B8]; "/daletqubuts", [0x05D3; 0x05BB];
"/daletqubutshebrew", [0x05D3; 0x05BB]; "/daletsegol", [0x05D3; 0x05B6];
"/daletsegolhebrew", [0x05D3; 0x05B6]; "/daletsheva", [0x05D3; 0x05B0];
"/daletshevahebrew", [0x05D3; 0x05B0]; "/dalettsere", [0x05D3; 0x05B5];
"/dalettserehebrew", [0x05D3; 0x05B5]; "/dalfinalarabic", [0xFEAA];
"/dammaarabic", [0x064F]; "/dammalowarabic", [0x064F]; "/dammatanaltonearabic",
[0x064C]; "/dammatanarabic", [0x064C]; "/danda", [0x0964]; "/dargahebrew",
[0x05A7]; "/dargalefthebrew", [0x05A7]; "/dasiapneumatacyrilliccmb", [0x0485];
"/dblGrave", [0xF6D3]; "/dblanglebracketleft", [0x300A];
"/dblanglebracketleftvertical", [0xFE3D]; "/dblanglebracketright", [0x300B];
"/dblanglebracketrightvertical", [0xFE3E]; "/dblarchinvertedbelowcmb", [0x032B];
"/dblarrowleft", [0x21D4]; "/dblarrowright", [0x21D2]; "/dbldanda", [0x0965];
"/dblgrave", [0xF6D6]; "/dblgravecmb", [0x030F]; "/dblintegral", [0x222C];
"/dbllowline", [0x2017]; "/dbllowlinecmb", [0x0333]; "/dbloverlinecmb",
[0x033F]; "/dblprimemod", [0x02BA]; "/dblverticalbar", [0x2016];
"/dblverticallineabovecmb", [0x030E]; "/dbopomofo", [0x3109]; "/dbsquare",
[0x33C8]; "/dcaron", [0x010F]; "/dcedilla", [0x1E11]; "/dcircle", [0x24D3];
"/dcircumflexbelow", [0x1E13]; "/dcroat", [0x0111]; "/ddabengali", [0x09A1];
"/ddadeva", [0x0921]; "/ddagujarati", [0x0AA1]; "/ddagurmukhi", [0x0A21];
"/ddalarabic", [0x0688]; "/ddalfinalarabic", [0xFB89]; "/dddhadeva", [0x095C];
"/ddhabengali", [0x09A2]; "/ddhadeva", [0x0922]; "/ddhagujarati", [0x0AA2];
"/ddhagurmukhi", [0x0A22]; "/ddotaccent", [0x1E0B]; "/ddotbelow", [0x1E0D];
"/decimalseparatorarabic", [0x066B]; "/decimalseparatorpersian", [0x066B];
"/decyrillic", [0x0434]; "/degree", [0x00B0]; "/dehihebrew", [0x05AD];
"/dehiragana", [0x3067]; "/deicoptic", [0x03EF]; "/dekatakana", [0x30C7];
"/deleteleft", [0x232B]; "/deleteright", [0x2326]; "/delta", [0x03B4];
"/deltaturned", [0x018D]; "/denominatorminusonenumeratorbengali", [0x09F8];
"/dezh", [0x02A4]; "/dhabengali", [0x09A7]; "/dhadeva", [0x0927];
"/dhagujarati", [0x0AA7]; "/dhagurmukhi", [0x0A27]; "/dhook", [0x0257];
"/dialytikatonos", [0x0385]; "/dialytikatonoscmb", [0x0344]; "/diamond",
[0x2666]; "/diamondsuitwhite", [0x2662]; "/dieresis", [0x00A8];
"/dieresisacute", [0xF6D7]; "/dieresisbelowcmb", [0x0324]; "/dieresiscmb",
[0x0308]; "/dieresisgrave", [0xF6D8]; "/dieresistonos", [0x0385]; "/dihiragana",
[0x3062]; "/dikatakana", [0x30C2]; "/dittomark", [0x3003]; "/divide", [0x00F7];
"/divides", [0x2223]; "/divisionslash", [0x2215]; "/djecyrillic", [0x0452];
"/dkshade", [0x2593]; "/dlinebelow", [0x1E0F]; "/dlsquare", [0x3397];
"/dmacron", [0x0111]; "/dmonospace", [0xFF44]; "/dnblock", [0x2584];
"/dochadathai", [0x0E0E]; "/dodekthai", [0x0E14]; "/dohiragana", [0x3069];
"/dokatakana", [0x30C9]; "/dollar", [0x0024]; "/dollarinferior", [0xF6E3];
"/dollarmonospace", [0xFF04]; "/dollaroldstyle", [0xF724]; "/dollarsmall",
[0xFE69]; "/dollarsuperior", [0xF6E4]; "/dong", [0x20AB]; "/dorusquare",
[0x3326]; "/dotaccent", [0x02D9]; "/dotaccentcmb", [0x0307]; "/dotbelowcmb",
[0x0323]; "/dotbelowcomb", [0x0323]; "/dotkatakana", [0x30FB]; "/dotlessi",
[0x0131]; "/dotlessj", [0xF6BE]; "/dotlessjstrokehook", [0x0284]; "/dotmath",
[0x22C5]; "/dottedcircle", [0x25CC]; "/doubleyodpatah", [0xFB1F];
"/doubleyodpatahhebrew", [0xFB1F]; "/downtackbelowcmb", [0x031E];
"/downtackmod", [0x02D5]; "/dparen", [0x249F]; "/dsuperior", [0xF6EB]; "/dtail",
[0x0256]; "/dtopbar", [0x018C]; "/duhiragana", [0x3065]; "/dukatakana",
[0x30C5]; "/dz", [0x01F3]; "/dzaltone", [0x02A3]; "/dzcaron", [0x01C6];
"/dzcurl", [0x02A5]; "/dzeabkhasiancyrillic", [0x04E1]; "/dzecyrillic",
[0x0455]; "/dzhecyrillic", [0x045F]; "/e", [0x0065]; "/eacute", [0x00E9];
"/earth", [0x2641]; "/ebengali", [0x098F]; "/ebopomofo", [0x311C]; "/ebreve",
[0x0115]; "/ecandradeva", [0x090D]; "/ecandragujarati", [0x0A8D];
"/ecandravowelsigndeva", [0x0945]; "/ecandravowelsigngujarati", [0x0AC5];
"/ecaron", [0x011B]; "/ecedillabreve", [0x1E1D]; "/echarmenian", [0x0565];
"/echyiwnarmenian", [0x0587]; "/ecircle", [0x24D4]; "/ecircumflex", [0x00EA];
"/ecircumflexacute", [0x1EBF]; "/ecircumflexbelow", [0x1E19];
"/ecircumflexdotbelow", [0x1EC7]; "/ecircumflexgrave", [0x1EC1];
"/ecircumflexhookabove", [0x1EC3]; "/ecircumflextilde", [0x1EC5]; "/ecyrillic",
[0x0454]; "/edblgrave", [0x0205]; "/edeva", [0x090F]; "/edieresis", [0x00EB];
"/edot", [0x0117]; "/edotaccent", [0x0117]; "/edotbelow", [0x1EB9];
"/eegurmukhi", [0x0A0F]; "/eematragurmukhi", [0x0A47]; "/efcyrillic", [0x0444];
"/egrave", [0x00E8]; "/egujarati", [0x0A8F]; "/eharmenian", [0x0567];
"/ehbopomofo", [0x311D]; "/ehiragana", [0x3048]; "/ehookabove", [0x1EBB];
"/eibopomofo", [0x311F]; "/eight", [0x0038]; "/eightarabic", [0x0668];
"/eightbengali", [0x09EE]; "/eightcircle", [0x2467];
"/eightcircleinversesansserif", [0x2791]; "/eightdeva", [0x096E];
"/eighteencircle", [0x2471]; "/eighteenparen", [0x2485]; "/eighteenperiod",
[0x2499]; "/eightgujarati", [0x0AEE]; "/eightgurmukhi", [0x0A6E];
"/eighthackarabic", [0x0668]; "/eighthangzhou", [0x3028]; "/eighthnotebeamed",
[0x266B]; "/eightideographicparen", [0x3227]; "/eightinferior", [0x2088];
"/eightmonospace", [0xFF18]; "/eightoldstyle", [0xF738]; "/eightparen",
[0x247B]; "/eightperiod", [0x248F]; "/eightpersian", [0x06F8]; "/eightroman",
[0x2177]; "/eightsuperior", [0x2078]; "/eightthai", [0x0E58]; "/einvertedbreve",
[0x0207]; "/eiotifiedcyrillic", [0x0465]; "/ekatakana", [0x30A8];
"/ekatakanahalfwidth", [0xFF74]; "/ekonkargurmukhi", [0x0A74]; "/ekorean",
[0x3154]; "/elcyrillic", [0x043B]; "/element", [0x2208]; "/elevencircle",
[0x246A]; "/elevenparen", [0x247E]; "/elevenperiod", [0x2492]; "/elevenroman",
[0x217A]; "/ellipsis", [0x2026]; "/ellipsisvertical", [0x22EE]; "/emacron",
[0x0113]; "/emacronacute", [0x1E17]; "/emacrongrave", [0x1E15]; "/emcyrillic",
[0x043C]; "/emdash", [0x2014]; "/emdashvertical", [0xFE31]; "/emonospace",
[0xFF45]; "/emphasismarkarmenian", [0x055B]; "/emptyset", [0x2205];
"/enbopomofo", [0x3123]; "/encyrillic", [0x043D]; "/endash", [0x2013];
"/endashvertical", [0xFE32]; "/endescendercyrillic", [0x04A3]; "/eng", [0x014B];
"/engbopomofo", [0x3125]; "/enghecyrillic", [0x04A5]; "/enhookcyrillic",
[0x04C8]; "/enspace", [0x2002]; "/eogonek", [0x0119]; "/eokorean", [0x3153];
"/eopen", [0x025B]; "/eopenclosed", [0x029A]; "/eopenreversed", [0x025C];
"/eopenreversedclosed", [0x025E]; "/eopenreversedhook", [0x025D]; "/eparen",
[0x24A0]; "/epsilon", [0x03B5]; "/epsilontonos", [0x03AD]; "/equal", [0x003D];
"/equalmonospace", [0xFF1D]; "/equalsmall", [0xFE66]; "/equalsuperior",
[0x207C]; "/equivalence", [0x2261]; "/erbopomofo", [0x3126]; "/ercyrillic",
[0x0440]; "/ereversed", [0x0258]; "/ereversedcyrillic", [0x044D]; "/escyrillic",
[0x0441]; "/esdescendercyrillic", [0x04AB]; "/esh", [0x0283]; "/eshcurl",
[0x0286]; "/eshortdeva", [0x090E]; "/eshortvowelsigndeva", [0x0946];
"/eshreversedloop", [0x01AA]; "/eshsquatreversed", [0x0285]; "/esmallhiragana",
[0x3047]; "/esmallkatakana", [0x30A7]; "/esmallkatakanahalfwidth", [0xFF6A];
"/estimated", [0x212E]; "/esuperior", [0xF6EC]; "/eta", [0x03B7]; "/etarmenian",
[0x0568]; "/etatonos", [0x03AE]; "/eth", [0x00F0]; "/etilde", [0x1EBD];
"/etildebelow", [0x1E1B]; "/etnahtafoukhhebrew", [0x0591];
"/etnahtafoukhlefthebrew", [0x0591]; "/etnahtahebrew", [0x0591];
"/etnahtalefthebrew", [0x0591]; "/eturned", [0x01DD]; "/eukorean", [0x3161];
"/euro", [0x20AC]; "/evowelsignbengali", [0x09C7]; "/evowelsigndeva", [0x0947];
"/evowelsigngujarati", [0x0AC7]; "/exclam", [0x0021]; "/exclamarmenian",
[0x055C]; "/exclamdbl", [0x203C]; "/exclamdown", [0x00A1]; "/exclamdownsmall",
[0xF7A1]; "/exclammonospace", [0xFF01]; "/exclamsmall", [0xF721];
"/existential", [0x2203]; "/ezh", [0x0292]; "/ezhcaron", [0x01EF]; "/ezhcurl",
[0x0293]; "/ezhreversed", [0x01B9]; "/ezhtail", [0x01BA]; "/f", [0x0066];
"/fadeva", [0x095E]; "/fagurmukhi", [0x0A5E]; "/fahrenheit", [0x2109];
"/fathaarabic", [0x064E]; "/fathalowarabic", [0x064E]; "/fathatanarabic",
[0x064B]; "/fbopomofo", [0x3108]; "/fcircle", [0x24D5]; "/fdotaccent", [0x1E1F];
"/feharabic", [0x0641]; "/feharmenian", [0x0586]; "/fehfinalarabic", [0xFED2];
"/fehinitialarabic", [0xFED3]; "/fehmedialarabic", [0xFED4]; "/feicoptic",
[0x03E5]; "/female", [0x2640]; "/ff", [0xFB00]; "/ffi", [0xFB03]; "/ffl",
[0xFB04]; "/fi", [0xFB01]; "/fifteencircle", [0x246E]; "/fifteenparen",
[0x2482]; "/fifteenperiod", [0x2496]; "/figuredash", [0x2012]; "/filledbox",
[0x25A0]; "/filledrect", [0x25AC]; "/finalkaf", [0x05DA]; "/finalkafdagesh",
[0xFB3A]; "/finalkafdageshhebrew", [0xFB3A]; "/finalkafhebrew", [0x05DA];
"/finalkafqamats", [0x05DA; 0x05B8]; "/finalkafqamatshebrew", [0x05DA; 0x05B8];
"/finalkafsheva", [0x05DA; 0x05B0]; "/finalkafshevahebrew", [0x05DA; 0x05B0];
"/finalmem", [0x05DD]; "/finalmemhebrew", [0x05DD]; "/finalnun", [0x05DF];
"/finalnunhebrew", [0x05DF]; "/finalpe", [0x05E3]; "/finalpehebrew", [0x05E3];
"/finaltsadi", [0x05E5]; "/finaltsadihebrew", [0x05E5]; "/firsttonechinese",
[0x02C9]; "/fisheye", [0x25C9]; "/fitacyrillic", [0x0473]; "/five", [0x0035];
"/fivearabic", [0x0665]; "/fivebengali", [0x09EB]; "/fivecircle", [0x2464];
"/fivecircleinversesansserif", [0x278E]; "/fivedeva", [0x096B]; "/fiveeighths",
[0x215D]; "/fivegujarati", [0x0AEB]; "/fivegurmukhi", [0x0A6B];
"/fivehackarabic", [0x0665]; "/fivehangzhou", [0x3025]; "/fiveideographicparen",
[0x3224]; "/fiveinferior", [0x2085]; "/fivemonospace", [0xFF15];
"/fiveoldstyle", [0xF735]; "/fiveparen", [0x2478]; "/fiveperiod", [0x248C];
"/fivepersian", [0x06F5]; "/fiveroman", [0x2174]; "/fivesuperior", [0x2075];
"/fivethai", [0x0E55]; "/fl", [0xFB02]; "/florin", [0x0192]; "/fmonospace",
[0xFF46]; "/fmsquare", [0x3399]; "/fofanthai", [0x0E1F]; "/fofathai", [0x0E1D];
"/fongmanthai", [0x0E4F]; "/forall", [0x2200]; "/four", [0x0034]; "/fourarabic",
[0x0664]; "/fourbengali", [0x09EA]; "/fourcircle", [0x2463];
"/fourcircleinversesansserif", [0x278D]; "/fourdeva", [0x096A]; "/fourgujarati",
[0x0AEA]; "/fourgurmukhi", [0x0A6A]; "/fourhackarabic", [0x0664];
"/fourhangzhou", [0x3024]; "/fourideographicparen", [0x3223]; "/fourinferior",
[0x2084]; "/fourmonospace", [0xFF14]; "/fournumeratorbengali", [0x09F7];
"/fouroldstyle", [0xF734]; "/fourparen", [0x2477]; "/fourperiod", [0x248B];
"/fourpersian", [0x06F4]; "/fourroman", [0x2173]; "/foursuperior", [0x2074];
"/fourteencircle", [0x246D]; "/fourteenparen", [0x2481]; "/fourteenperiod",
[0x2495]; "/fourthai", [0x0E54]; "/fourthtonechinese", [0x02CB]; "/fparen",
[0x24A1]; "/fraction", [0x2044]; "/franc", [0x20A3]; "/g", [0x0067];
"/gabengali", [0x0997]; "/gacute", [0x01F5]; "/gadeva", [0x0917]; "/gafarabic",
[0x06AF]; "/gaffinalarabic", [0xFB93]; "/gafinitialarabic", [0xFB94];
"/gafmedialarabic", [0xFB95]; "/gagujarati", [0x0A97]; "/gagurmukhi", [0x0A17];
"/gahiragana", [0x304C]; "/gakatakana", [0x30AC]; "/gamma", [0x03B3];
"/gammalatinsmall", [0x0263]; "/gammasuperior", [0x02E0]; "/gangiacoptic",
[0x03EB]; "/gbopomofo", [0x310D]; "/gbreve", [0x011F]; "/gcaron", [0x01E7];
"/gcedilla", [0x0123]; "/gcircle", [0x24D6]; "/gcircumflex", [0x011D];
"/gcommaaccent", [0x0123]; "/gdot", [0x0121]; "/gdotaccent", [0x0121];
"/gecyrillic", [0x0433]; "/gehiragana", [0x3052]; "/gekatakana", [0x30B2];
"/geometricallyequal", [0x2251]; "/gereshaccenthebrew", [0x059C];
"/gereshhebrew", [0x05F3]; "/gereshmuqdamhebrew", [0x059D]; "/germandbls",
[0x00DF]; "/gershayimaccenthebrew", [0x059E]; "/gershayimhebrew", [0x05F4];
"/getamark", [0x3013]; "/ghabengali", [0x0998]; "/ghadarmenian", [0x0572];
"/ghadeva", [0x0918]; "/ghagujarati", [0x0A98]; "/ghagurmukhi", [0x0A18];
"/ghainarabic", [0x063A]; "/ghainfinalarabic", [0xFECE]; "/ghaininitialarabic",
[0xFECF]; "/ghainmedialarabic", [0xFED0]; "/ghemiddlehookcyrillic", [0x0495];
"/ghestrokecyrillic", [0x0493]; "/gheupturncyrillic", [0x0491]; "/ghhadeva",
[0x095A]; "/ghhagurmukhi", [0x0A5A]; "/ghook", [0x0260]; "/ghzsquare", [0x3393];
"/gihiragana", [0x304E]; "/gikatakana", [0x30AE]; "/gimarmenian", [0x0563];
"/gimel", [0x05D2]; "/gimeldagesh", [0xFB32]; "/gimeldageshhebrew", [0xFB32];
"/gimelhebrew", [0x05D2]; "/gjecyrillic", [0x0453]; "/glottalinvertedstroke",
[0x01BE]; "/glottalstop", [0x0294]; "/glottalstopinverted", [0x0296];
"/glottalstopmod", [0x02C0]; "/glottalstopreversed", [0x0295];
"/glottalstopreversedmod", [0x02C1]; "/glottalstopreversedsuperior", [0x02E4];
"/glottalstopstroke", [0x02A1]; "/glottalstopstrokereversed", [0x02A2];
"/gmacron", [0x1E21]; "/gmonospace", [0xFF47]; "/gohiragana", [0x3054];
"/gokatakana", [0x30B4]; "/gparen", [0x24A2]; "/gpasquare", [0x33AC];
"/gradient", [0x2207]; "/grave", [0x0060]; "/gravebelowcmb", [0x0316];
"/gravecmb", [0x0300]; "/gravecomb", [0x0300]; "/gravedeva", [0x0953];
"/gravelowmod", [0x02CE]; "/gravemonospace", [0xFF40]; "/gravetonecmb",
[0x0340]; "/greater", [0x003E]; "/greaterequal", [0x2265];
"/greaterequalorless", [0x22DB]; "/greatermonospace", [0xFF1E];
"/greaterorequivalent", [0x2273]; "/greaterorless", [0x2277];
"/greateroverequal", [0x2267]; "/greatersmall", [0xFE65]; "/gscript", [0x0261];
"/gstroke", [0x01E5]; "/guhiragana", [0x3050]; "/guillemotleft", [0x00AB];
"/guillemotright", [0x00BB]; "/guilsinglleft", [0x2039]; "/guilsinglright",
[0x203A]; "/gukatakana", [0x30B0]; "/guramusquare", [0x3318]; "/gysquare",
[0x33C9]; "/h", [0x0068]; "/haabkhasiancyrillic", [0x04A9]; "/haaltonearabic",
[0x06C1]; "/habengali", [0x09B9]; "/hadescendercyrillic", [0x04B3]; "/hadeva",
[0x0939]; "/hagujarati", [0x0AB9]; "/hagurmukhi", [0x0A39]; "/haharabic",
[0x062D]; "/hahfinalarabic", [0xFEA2]; "/hahinitialarabic", [0xFEA3];
"/hahiragana", [0x306F]; "/hahmedialarabic", [0xFEA4]; "/haitusquare", [0x332A];
"/hakatakana", [0x30CF]; "/hakatakanahalfwidth", [0xFF8A]; "/halantgurmukhi",
[0x0A4D]; "/hamzaarabic", [0x0621]; "/hamzadammaarabic", [0x0621; 0x064F];
"/hamzadammatanarabic", [0x0621; 0x064C]; "/hamzafathaarabic", [0x0621; 0x064E];
"/hamzafathatanarabic", [0x0621; 0x064B]; "/hamzalowarabic", [0x0621];
"/hamzalowkasraarabic", [0x0621; 0x0650]; "/hamzalowkasratanarabic", [0x0621; 0x064D];
"/hamzasukunarabic", [0x0621; 0x0652]; "/hangulfiller", [0x3164];
"/hardsigncyrillic", [0x044A]; "/harpoonleftbarbup", [0x21BC];
"/harpoonrightbarbup", [0x21C0]; "/hasquare", [0x33CA]; "/hatafpatah", [0x05B2];
"/hatafpatah16", [0x05B2]; "/hatafpatah23", [0x05B2]; "/hatafpatah2f", [0x05B2];
"/hatafpatahhebrew", [0x05B2]; "/hatafpatahnarrowhebrew", [0x05B2];
"/hatafpatahquarterhebrew", [0x05B2]; "/hatafpatahwidehebrew", [0x05B2];
"/hatafqamats", [0x05B3]; "/hatafqamats1b", [0x05B3]; "/hatafqamats28",
[0x05B3]; "/hatafqamats34", [0x05B3]; "/hatafqamatshebrew", [0x05B3];
"/hatafqamatsnarrowhebrew", [0x05B3]; "/hatafqamatsquarterhebrew", [0x05B3];
"/hatafqamatswidehebrew", [0x05B3]; "/hatafsegol", [0x05B1]; "/hatafsegol17",
[0x05B1]; "/hatafsegol24", [0x05B1]; "/hatafsegol30", [0x05B1];
"/hatafsegolhebrew", [0x05B1]; "/hatafsegolnarrowhebrew", [0x05B1];
"/hatafsegolquarterhebrew", [0x05B1]; "/hatafsegolwidehebrew", [0x05B1];
"/hbar", [0x0127]; "/hbopomofo", [0x310F]; "/hbrevebelow", [0x1E2B];
"/hcedilla", [0x1E29]; "/hcircle", [0x24D7]; "/hcircumflex", [0x0125];
"/hdieresis", [0x1E27]; "/hdotaccent", [0x1E23]; "/hdotbelow", [0x1E25]; "/he",
[0x05D4]; "/heart", [0x2665]; "/heartsuitblack", [0x2665]; "/heartsuitwhite",
[0x2661]; "/hedagesh", [0xFB34]; "/hedageshhebrew", [0xFB34];
"/hehaltonearabic", [0x06C1]; "/heharabic", [0x0647]; "/hehebrew", [0x05D4];
"/hehfinalaltonearabic", [0xFBA7]; "/hehfinalalttwoarabic", [0xFEEA];
"/hehfinalarabic", [0xFEEA]; "/hehhamzaabovefinalarabic", [0xFBA5];
"/hehhamzaaboveisolatedarabic", [0xFBA4]; "/hehinitialaltonearabic", [0xFBA8];
"/hehinitialarabic", [0xFEEB]; "/hehiragana", [0x3078];
"/hehmedialaltonearabic", [0xFBA9]; "/hehmedialarabic", [0xFEEC];
"/heiseierasquare", [0x337B]; "/hekatakana", [0x30D8]; "/hekatakanahalfwidth",
[0xFF8D]; "/hekutaarusquare", [0x3336]; "/henghook", [0x0267]; "/herutusquare",
[0x3339]; "/het", [0x05D7]; "/hethebrew", [0x05D7]; "/hhook", [0x0266];
"/hhooksuperior", [0x02B1]; "/hieuhacirclekorean", [0x327B];
"/hieuhaparenkorean", [0x321B]; "/hieuhcirclekorean", [0x326D]; "/hieuhkorean",
[0x314E]; "/hieuhparenkorean", [0x320D]; "/hihiragana", [0x3072]; "/hikatakana",
[0x30D2]; "/hikatakanahalfwidth", [0xFF8B]; "/hiriq", [0x05B4]; "/hiriq14",
[0x05B4]; "/hiriq21", [0x05B4]; "/hiriq2d", [0x05B4]; "/hiriqhebrew", [0x05B4];
"/hiriqnarrowhebrew", [0x05B4]; "/hiriqquarterhebrew", [0x05B4];
"/hiriqwidehebrew", [0x05B4]; "/hlinebelow", [0x1E96]; "/hmonospace", [0xFF48];
"/hoarmenian", [0x0570]; "/hohipthai", [0x0E2B]; "/hohiragana", [0x307B];
"/hokatakana", [0x30DB]; "/hokatakanahalfwidth", [0xFF8E]; "/holam", [0x05B9];
"/holam19", [0x05B9]; "/holam26", [0x05B9]; "/holam32", [0x05B9];
"/holamhebrew", [0x05B9]; "/holamnarrowhebrew", [0x05B9]; "/holamquarterhebrew",
[0x05B9]; "/holamwidehebrew", [0x05B9]; "/honokhukthai", [0x0E2E];
"/hookabovecomb", [0x0309]; "/hookcmb", [0x0309]; "/hookpalatalizedbelowcmb",
[0x0321]; "/hookretroflexbelowcmb", [0x0322]; "/hoonsquare", [0x3342];
"/horicoptic", [0x03E9]; "/horizontalbar", [0x2015]; "/horncmb", [0x031B];
"/hotsprings", [0x2668]; "/house", [0x2302]; "/hparen", [0x24A3]; "/hsuperior",
[0x02B0]; "/hturned", [0x0265]; "/huhiragana", [0x3075]; "/huiitosquare",
[0x3333]; "/hukatakana", [0x30D5]; "/hukatakanahalfwidth", [0xFF8C];
"/hungarumlaut", [0x02DD]; "/hungarumlautcmb", [0x030B]; "/hv", [0x0195];
"/hyphen", [0x002D]; "/hypheninferior", [0xF6E5]; "/hyphenmonospace", [0xFF0D];
"/hyphensmall", [0xFE63]; "/hyphensuperior", [0xF6E6]; "/hyphentwo", [0x2010];
"/i", [0x0069]; "/iacute", [0x00ED]; "/iacyrillic", [0x044F]; "/ibengali",
[0x0987]; "/ibopomofo", [0x3127]; "/ibreve", [0x012D]; "/icaron", [0x01D0];
"/icircle", [0x24D8]; "/icircumflex", [0x00EE]; "/icyrillic", [0x0456];
"/idblgrave", [0x0209]; "/ideographearthcircle", [0x328F];
"/ideographfirecircle", [0x328B]; "/ideographicallianceparen", [0x323F];
"/ideographiccallparen", [0x323A]; "/ideographiccentrecircle", [0x32A5];
"/ideographicclose", [0x3006]; "/ideographiccomma", [0x3001];
"/ideographiccommaleft", [0xFF64]; "/ideographiccongratulationparen", [0x3237];
"/ideographiccorrectcircle", [0x32A3]; "/ideographicearthparen", [0x322F];
"/ideographicenterpriseparen", [0x323D]; "/ideographicexcellentcircle",
[0x329D]; "/ideographicfestivalparen", [0x3240]; "/ideographicfinancialcircle",
[0x3296]; "/ideographicfinancialparen", [0x3236]; "/ideographicfireparen",
[0x322B]; "/ideographichaveparen", [0x3232]; "/ideographichighcircle", [0x32A4];
"/ideographiciterationmark", [0x3005]; "/ideographiclaborcircle", [0x3298];
"/ideographiclaborparen", [0x3238]; "/ideographicleftcircle", [0x32A7];
"/ideographiclowcircle", [0x32A6]; "/ideographicmedicinecircle", [0x32A9];
"/ideographicmetalparen", [0x322E]; "/ideographicmoonparen", [0x322A];
"/ideographicnameparen", [0x3234]; "/ideographicperiod", [0x3002];
"/ideographicprintcircle", [0x329E]; "/ideographicreachparen", [0x3243];
"/ideographicrepresentparen", [0x3239]; "/ideographicresourceparen", [0x323E];
"/ideographicrightcircle", [0x32A8]; "/ideographicsecretcircle", [0x3299];
"/ideographicselfparen", [0x3242]; "/ideographicsocietyparen", [0x3233];
"/ideographicspace", [0x3000]; "/ideographicspecialparen", [0x3235];
"/ideographicstockparen", [0x3231]; "/ideographicstudyparen", [0x323B];
"/ideographicsunparen", [0x3230]; "/ideographicsuperviseparen", [0x323C];
"/ideographicwaterparen", [0x322C]; "/ideographicwoodparen", [0x322D];
"/ideographiczero", [0x3007]; "/ideographmetalcircle", [0x328E];
"/ideographmooncircle", [0x328A]; "/ideographnamecircle", [0x3294];
"/ideographsuncircle", [0x3290]; "/ideographwatercircle", [0x328C];
"/ideographwoodcircle", [0x328D]; "/ideva", [0x0907]; "/idieresis", [0x00EF];
"/idieresisacute", [0x1E2F]; "/idieresiscyrillic", [0x04E5]; "/idotbelow",
[0x1ECB]; "/iebrevecyrillic", [0x04D7]; "/iecyrillic", [0x0435];
"/ieungacirclekorean", [0x3275]; "/ieungaparenkorean", [0x3215];
"/ieungcirclekorean", [0x3267]; "/ieungkorean", [0x3147]; "/ieungparenkorean",
[0x3207]; "/igrave", [0x00EC]; "/igujarati", [0x0A87]; "/igurmukhi", [0x0A07];
"/ihiragana", [0x3044]; "/ihookabove", [0x1EC9]; "/iibengali", [0x0988];
"/iicyrillic", [0x0438]; "/iideva", [0x0908]; "/iigujarati", [0x0A88];
"/iigurmukhi", [0x0A08]; "/iimatragurmukhi", [0x0A40]; "/iinvertedbreve",
[0x020B]; "/iishortcyrillic", [0x0439]; "/iivowelsignbengali", [0x09C0];
"/iivowelsigndeva", [0x0940]; "/iivowelsigngujarati", [0x0AC0]; "/ij", [0x0133];
"/ikatakana", [0x30A4]; "/ikatakanahalfwidth", [0xFF72]; "/ikorean", [0x3163];
"/ilde", [0x02DC]; "/iluyhebrew", [0x05AC]; "/imacron", [0x012B];
"/imacroncyrillic", [0x04E3]; "/imageorapproximatelyequal", [0x2253];
"/imatragurmukhi", [0x0A3F]; "/imonospace", [0xFF49]; "/increment", [0x2206];
"/infinity", [0x221E]; "/iniarmenian", [0x056B]; "/integral", [0x222B];
"/integralbottom", [0x2321]; "/integralbt", [0x2321]; "/integralex", [0xF8F5];
"/integraltop", [0x2320]; "/integraltp", [0x2320]; "/intersection", [0x2229];
"/intisquare", [0x3305]; "/invbullet", [0x25D8]; "/invcircle", [0x25D9];
"/invsmileface", [0x263B]; "/iocyrillic", [0x0451]; "/iogonek", [0x012F];
"/iota", [0x03B9]; "/iotadieresis", [0x03CA]; "/iotadieresistonos", [0x0390];
"/iotalatin", [0x0269]; "/iotatonos", [0x03AF]; "/iparen", [0x24A4];
"/irigurmukhi", [0x0A72]; "/ismallhiragana", [0x3043]; "/ismallkatakana",
[0x30A3]; "/ismallkatakanahalfwidth", [0xFF68]; "/issharbengali", [0x09FA];
"/istroke", [0x0268]; "/isuperior", [0xF6ED]; "/iterationhiragana", [0x309D];
"/iterationkatakana", [0x30FD]; "/itilde", [0x0129]; "/itildebelow", [0x1E2D];
"/iubopomofo", [0x3129]; "/iucyrillic", [0x044E]; "/ivowelsignbengali",
[0x09BF]; "/ivowelsigndeva", [0x093F]; "/ivowelsigngujarati", [0x0ABF];
"/izhitsacyrillic", [0x0475]; "/izhitsadblgravecyrillic", [0x0477]; "/j",
[0x006A]; "/jaarmenian", [0x0571]; "/jabengali", [0x099C]; "/jadeva", [0x091C];
"/jagujarati", [0x0A9C]; "/jagurmukhi", [0x0A1C]; "/jbopomofo", [0x3110];
"/jcaron", [0x01F0]; "/jcircle", [0x24D9]; "/jcircumflex", [0x0135];
"/jcrossedtail", [0x029D]; "/jdotlessstroke", [0x025F]; "/jecyrillic", [0x0458];
"/jeemarabic", [0x062C]; "/jeemfinalarabic", [0xFE9E]; "/jeeminitialarabic",
[0xFE9F]; "/jeemmedialarabic", [0xFEA0]; "/jeharabic", [0x0698];
"/jehfinalarabic", [0xFB8B]; "/jhabengali", [0x099D]; "/jhadeva", [0x091D];
"/jhagujarati", [0x0A9D]; "/jhagurmukhi", [0x0A1D]; "/jheharmenian", [0x057B];
"/jis", [0x3004]; "/jmonospace", [0xFF4A]; "/jparen", [0x24A5]; "/jsuperior",
[0x02B2]; "/k", [0x006B]; "/kabashkircyrillic", [0x04A1]; "/kabengali",
[0x0995]; "/kacute", [0x1E31]; "/kacyrillic", [0x043A]; "/kadescendercyrillic",
[0x049B]; "/kadeva", [0x0915]; "/kaf", [0x05DB]; "/kafarabic", [0x0643];
"/kafdagesh", [0xFB3B]; "/kafdageshhebrew", [0xFB3B]; "/kaffinalarabic",
[0xFEDA]; "/kafhebrew", [0x05DB]; "/kafinitialarabic", [0xFEDB];
"/kafmedialarabic", [0xFEDC]; "/kafrafehebrew", [0xFB4D]; "/kagujarati",
[0x0A95]; "/kagurmukhi", [0x0A15]; "/kahiragana", [0x304B]; "/kahookcyrillic",
[0x04C4]; "/kakatakana", [0x30AB]; "/kakatakanahalfwidth", [0xFF76]; "/kappa",
[0x03BA]; "/kappasymbolgreek", [0x03F0]; "/kapyeounmieumkorean", [0x3171];
"/kapyeounphieuphkorean", [0x3184]; "/kapyeounpieupkorean", [0x3178];
"/kapyeounssangpieupkorean", [0x3179]; "/karoriisquare", [0x330D];
"/kashidaautoarabic", [0x0640]; "/kashidaautonosidebearingarabic", [0x0640];
"/kasmallkatakana", [0x30F5]; "/kasquare", [0x3384]; "/kasraarabic", [0x0650];
"/kasratanarabic", [0x064D]; "/kastrokecyrillic", [0x049F];
"/katahiraprolongmarkhalfwidth", [0xFF70]; "/kaverticalstrokecyrillic",
[0x049D]; "/kbopomofo", [0x310E]; "/kcalsquare", [0x3389]; "/kcaron", [0x01E9];
"/kcedilla", [0x0137]; "/kcircle", [0x24DA]; "/kcommaaccent", [0x0137];
"/kdotbelow", [0x1E33]; "/keharmenian", [0x0584]; "/kehiragana", [0x3051];
"/kekatakana", [0x30B1]; "/kekatakanahalfwidth", [0xFF79]; "/kenarmenian",
[0x056F]; "/kesmallkatakana", [0x30F6]; "/kgreenlandic", [0x0138];
"/khabengali", [0x0996]; "/khacyrillic", [0x0445]; "/khadeva", [0x0916];
"/khagujarati", [0x0A96]; "/khagurmukhi", [0x0A16]; "/khaharabic", [0x062E];
"/khahfinalarabic", [0xFEA6]; "/khahinitialarabic", [0xFEA7];
"/khahmedialarabic", [0xFEA8]; "/kheicoptic", [0x03E7]; "/khhadeva", [0x0959];
"/khhagurmukhi", [0x0A59]; "/khieukhacirclekorean", [0x3278];
"/khieukhaparenkorean", [0x3218]; "/khieukhcirclekorean", [0x326A];
"/khieukhkorean", [0x314B]; "/khieukhparenkorean", [0x320A]; "/khokhaithai",
[0x0E02]; "/khokhonthai", [0x0E05]; "/khokhuatthai", [0x0E03]; "/khokhwaithai",
[0x0E04]; "/khomutthai", [0x0E5B]; "/khook", [0x0199]; "/khorakhangthai",
[0x0E06]; "/khzsquare", [0x3391]; "/kihiragana", [0x304D]; "/kikatakana",
[0x30AD]; "/kikatakanahalfwidth", [0xFF77]; "/kiroguramusquare", [0x3315];
"/kiromeetorusquare", [0x3316]; "/kirosquare", [0x3314]; "/kiyeokacirclekorean",
[0x326E]; "/kiyeokaparenkorean", [0x320E]; "/kiyeokcirclekorean", [0x3260];
"/kiyeokkorean", [0x3131]; "/kiyeokparenkorean", [0x3200]; "/kiyeoksioskorean",
[0x3133]; "/kjecyrillic", [0x045C]; "/klinebelow", [0x1E35]; "/klsquare",
[0x3398]; "/kmcubedsquare", [0x33A6]; "/kmonospace", [0xFF4B];
"/kmsquaredsquare", [0x33A2]; "/kohiragana", [0x3053]; "/kohmsquare", [0x33C0];
"/kokaithai", [0x0E01]; "/kokatakana", [0x30B3]; "/kokatakanahalfwidth",
[0xFF7A]; "/kooposquare", [0x331E]; "/koppacyrillic", [0x0481];
"/koreanstandardsymbol", [0x327F]; "/koroniscmb", [0x0343]; "/kparen", [0x24A6];
"/kpasquare", [0x33AA]; "/ksicyrillic", [0x046F]; "/ktsquare", [0x33CF];
"/kturned", [0x029E]; "/kuhiragana", [0x304F]; "/kukatakana", [0x30AF];
"/kukatakanahalfwidth", [0xFF78]; "/kvsquare", [0x33B8]; "/kwsquare", [0x33BE];
"/l", [0x006C]; "/labengali", [0x09B2]; "/lacute", [0x013A]; "/ladeva",
[0x0932]; "/lagujarati", [0x0AB2]; "/lagurmukhi", [0x0A32]; "/lakkhangyaothai",
[0x0E45]; "/lamaleffinalarabic", [0xFEFC]; "/lamalefhamzaabovefinalarabic",
[0xFEF8]; "/lamalefhamzaaboveisolatedarabic", [0xFEF7];
"/lamalefhamzabelowfinalarabic", [0xFEFA]; "/lamalefhamzabelowisolatedarabic",
[0xFEF9]; "/lamalefisolatedarabic", [0xFEFB]; "/lamalefmaddaabovefinalarabic",
[0xFEF6]; "/lamalefmaddaaboveisolatedarabic", [0xFEF5]; "/lamarabic", [0x0644];
"/lambda", [0x03BB]; "/lambdastroke", [0x019B]; "/lamed", [0x05DC];
"/lameddagesh", [0xFB3C]; "/lameddageshhebrew", [0xFB3C]; "/lamedhebrew",
[0x05DC]; "/lamedholam", [0x05DC; 0x05B9]; "/lamedholamdagesh", [0x05DC; 0x05B9;
0x05BC];
"/lamedholamdageshhebrew", [0x05DC; 0x05B9; 0x05BC]; "/lamedholamhebrew",
[0x05DC;
0x05B9]; "/lamfinalarabic", [0xFEDE]; "/lamhahinitialarabic", [0xFCCA];
"/laminitialarabic", [0xFEDF]; "/lamjeeminitialarabic", [0xFCC9];
"/lamkhahinitialarabic", [0xFCCB]; "/lamlamhehisolatedarabic", [0xFDF2];
"/lammedialarabic", [0xFEE0]; "/lammeemhahinitialarabic", [0xFD88];
"/lammeeminitialarabic", [0xFCCC]; "/lammeemjeeminitialarabic", [0xFEDF; 0xFEE4;
 0xFEA0]; "/lammeemkhahinitialarabic", [0xFEDF; 0xFEE4; 0xFEA8]; "/largecircle",
[0x25EF]; "/lbar", [0x019A]; "/lbelt", [0x026C]; "/lbopomofo", [0x310C];
"/lcaron", [0x013E]; "/lcedilla", [0x013C]; "/lcircle", [0x24DB];
"/lcircumflexbelow", [0x1E3D]; "/lcommaaccent", [0x013C]; "/ldot", [0x0140];
"/ldotaccent", [0x0140]; "/ldotbelow", [0x1E37]; "/ldotbelowmacron", [0x1E39];
"/leftangleabovecmb", [0x031A]; "/lefttackbelowcmb", [0x0318]; "/less",
[0x003C]; "/lessequal", [0x2264]; "/lessequalorgreater", [0x22DA];
"/lessmonospace", [0xFF1C]; "/lessorequivalent", [0x2272]; "/lessorgreater",
[0x2276]; "/lessoverequal", [0x2266]; "/lesssmall", [0xFE64]; "/lezh", [0x026E];
"/lfblock", [0x258C]; "/lhookretroflex", [0x026D]; "/lira", [0x20A4];
"/liwnarmenian", [0x056C]; "/lj", [0x01C9]; "/ljecyrillic", [0x0459]; "/ll",
[0xF6C0]; "/lladeva", [0x0933]; "/llagujarati", [0x0AB3]; "/llinebelow",
[0x1E3B]; "/llladeva", [0x0934]; "/llvocalicbengali", [0x09E1];
"/llvocalicdeva", [0x0961]; "/llvocalicvowelsignbengali", [0x09E3];
"/llvocalicvowelsigndeva", [0x0963]; "/lmiddletilde", [0x026B]; "/lmonospace",
[0xFF4C]; "/lmsquare", [0x33D0]; "/lochulathai", [0x0E2C]; "/logicaland",
[0x2227]; "/logicalnot", [0x00AC]; "/logicalnotreversed", [0x2310];
"/logicalor", [0x2228]; "/lolingthai", [0x0E25]; "/longs", [0x017F];
"/lowlinecenterline", [0xFE4E]; "/lowlinecmb", [0x0332]; "/lowlinedashed",
[0xFE4D]; "/lozenge", [0x25CA]; "/lparen", [0x24A7]; "/lslash", [0x0142];
"/lsquare", [0x2113]; "/lsuperior", [0xF6EE]; "/ltshade", [0x2591]; "/luthai",
[0x0E26]; "/lvocalicbengali", [0x098C]; "/lvocalicdeva", [0x090C];
"/lvocalicvowelsignbengali", [0x09E2]; "/lvocalicvowelsigndeva", [0x0962];
"/lxsquare", [0x33D3]; "/m", [0x006D]; "/mabengali", [0x09AE]; "/macron",
[0x00AF]; "/macronbelowcmb", [0x0331]; "/macroncmb", [0x0304]; "/macronlowmod",
[0x02CD]; "/macronmonospace", [0xFFE3]; "/macute", [0x1E3F]; "/madeva",
[0x092E]; "/magujarati", [0x0AAE]; "/magurmukhi", [0x0A2E]; "/mahapakhhebrew",
[0x05A4]; "/mahapakhlefthebrew", [0x05A4]; "/mahiragana", [0x307E];
"/maichattawalowleftthai", [0xF895]; "/maichattawalowrightthai", [0xF894];
"/maichattawathai", [0x0E4B]; "/maichattawaupperleftthai", [0xF893];
"/maieklowleftthai", [0xF88C]; "/maieklowrightthai", [0xF88B]; "/maiekthai",
[0x0E48]; "/maiekupperleftthai", [0xF88A]; "/maihanakatleftthai", [0xF884];
"/maihanakatthai", [0x0E31]; "/maitaikhuleftthai", [0xF889]; "/maitaikhuthai",
[0x0E47]; "/maitholowleftthai", [0xF88F]; "/maitholowrightthai", [0xF88E];
"/maithothai", [0x0E49]; "/maithoupperleftthai", [0xF88D]; "/maitrilowleftthai",
[0xF892]; "/maitrilowrightthai", [0xF891]; "/maitrithai", [0x0E4A];
"/maitriupperleftthai", [0xF890]; "/maiyamokthai", [0x0E46]; "/makatakana",
[0x30DE]; "/makatakanahalfwidth", [0xFF8F]; "/male", [0x2642]; "/mansyonsquare",
[0x3347]; "/maqafhebrew", [0x05BE]; "/mars", [0x2642]; "/masoracirclehebrew",
[0x05AF]; "/masquare", [0x3383]; "/mbopomofo", [0x3107]; "/mbsquare", [0x33D4];
"/mcircle", [0x24DC]; "/mcubedsquare", [0x33A5]; "/mdotaccent", [0x1E41];
"/mdotbelow", [0x1E43]; "/meemarabic", [0x0645]; "/meemfinalarabic", [0xFEE2];
"/meeminitialarabic", [0xFEE3]; "/meemmedialarabic", [0xFEE4];
"/meemmeeminitialarabic", [0xFCD1]; "/meemmeemisolatedarabic", [0xFC48];
"/meetorusquare", [0x334D]; "/mehiragana", [0x3081]; "/meizierasquare",
[0x337E]; "/mekatakana", [0x30E1]; "/mekatakanahalfwidth", [0xFF92]; "/mem",
[0x05DE]; "/memdagesh", [0xFB3E]; "/memdageshhebrew", [0xFB3E]; "/memhebrew",
[0x05DE]; "/menarmenian", [0x0574]; "/merkhahebrew", [0x05A5];
"/merkhakefulahebrew", [0x05A6]; "/merkhakefulalefthebrew", [0x05A6];
"/merkhalefthebrew", [0x05A5]; "/mhook", [0x0271]; "/mhzsquare", [0x3392];
"/middledotkatakanahalfwidth", [0xFF65]; "/middot", [0x00B7];
"/mieumacirclekorean", [0x3272]; "/mieumaparenkorean", [0x3212];
"/mieumcirclekorean", [0x3264]; "/mieumkorean", [0x3141]; "/mieumpansioskorean",
[0x3170]; "/mieumparenkorean", [0x3204]; "/mieumpieupkorean", [0x316E];
"/mieumsioskorean", [0x316F]; "/mihiragana", [0x307F]; "/mikatakana", [0x30DF];
"/mikatakanahalfwidth", [0xFF90]; "/minus", [0x2212]; "/minusbelowcmb",
[0x0320]; "/minuscircle", [0x2296]; "/minusmod", [0x02D7]; "/minusplus",
[0x2213]; "/minute", [0x2032]; "/miribaarusquare", [0x334A]; "/mirisquare",
[0x3349]; "/mlonglegturned", [0x0270]; "/mlsquare", [0x3396]; "/mmcubedsquare",
[0x33A3]; "/mmonospace", [0xFF4D]; "/mmsquaredsquare", [0x339F]; "/mohiragana",
[0x3082]; "/mohmsquare", [0x33C1]; "/mokatakana", [0x30E2];
"/mokatakanahalfwidth", [0xFF93]; "/molsquare", [0x33D6]; "/momathai", [0x0E21];
"/moverssquare", [0x33A7]; "/moverssquaredsquare", [0x33A8]; "/mparen",
[0x24A8]; "/mpasquare", [0x33AB]; "/mssquare", [0x33B3]; "/msuperior", [0xF6EF];
"/mturned", [0x026F]; "/mu", [0x00B5]; "/mu1", [0x00B5]; "/muasquare", [0x3382];
"/muchgreater", [0x226B]; "/muchless", [0x226A]; "/mufsquare", [0x338C];
"/mugreek", [0x03BC]; "/mugsquare", [0x338D]; "/muhiragana", [0x3080];
"/mukatakana", [0x30E0]; "/mukatakanahalfwidth", [0xFF91]; "/mulsquare",
[0x3395]; "/multiply", [0x00D7]; "/mumsquare", [0x339B]; "/munahhebrew",
[0x05A3]; "/munahlefthebrew", [0x05A3]; "/musicalnote", [0x266A];
"/musicalnotedbl", [0x266B]; "/musicflatsign", [0x266D]; "/musicsharpsign",
[0x266F]; "/mussquare", [0x33B2]; "/muvsquare", [0x33B6]; "/muwsquare",
[0x33BC]; "/mvmegasquare", [0x33B9]; "/mvsquare", [0x33B7]; "/mwmegasquare",
[0x33BF]; "/mwsquare", [0x33BD]; "/n", [0x006E]; "/nabengali", [0x09A8];
"/nabla", [0x2207]; "/nacute", [0x0144]; "/nadeva", [0x0928]; "/nagujarati",
[0x0AA8]; "/nagurmukhi", [0x0A28]; "/nahiragana", [0x306A]; "/nakatakana",
[0x30CA]; "/nakatakanahalfwidth", [0xFF85]; "/napostrophe", [0x0149];
"/nasquare", [0x3381]; "/nbopomofo", [0x310B]; "/nbspace", [0x00A0]; "/ncaron",
[0x0148]; "/ncedilla", [0x0146]; "/ncircle", [0x24DD]; "/ncircumflexbelow",
[0x1E4B]; "/ncommaaccent", [0x0146]; "/ndotaccent", [0x1E45]; "/ndotbelow",
[0x1E47]; "/nehiragana", [0x306D]; "/nekatakana", [0x30CD];
"/nekatakanahalfwidth", [0xFF88]; "/newsheqelsign", [0x20AA]; "/nfsquare",
[0x338B]; "/ngabengali", [0x0999]; "/ngadeva", [0x0919]; "/ngagujarati",
[0x0A99]; "/ngagurmukhi", [0x0A19]; "/ngonguthai", [0x0E07]; "/nhiragana",
[0x3093]; "/nhookleft", [0x0272]; "/nhookretroflex", [0x0273];
"/nieunacirclekorean", [0x326F]; "/nieunaparenkorean", [0x320F];
"/nieuncieuckorean", [0x3135]; "/nieuncirclekorean", [0x3261];
"/nieunhieuhkorean", [0x3136]; "/nieunkorean", [0x3134]; "/nieunpansioskorean",
[0x3168]; "/nieunparenkorean", [0x3201]; "/nieunsioskorean", [0x3167];
"/nieuntikeutkorean", [0x3166]; "/nihiragana", [0x306B]; "/nikatakana",
[0x30CB]; "/nikatakanahalfwidth", [0xFF86]; "/nikhahitleftthai", [0xF899];
"/nikhahitthai", [0x0E4D]; "/nine", [0x0039]; "/ninearabic", [0x0669];
"/ninebengali", [0x09EF]; "/ninecircle", [0x2468];
"/ninecircleinversesansserif", [0x2792]; "/ninedeva", [0x096F]; "/ninegujarati",
[0x0AEF]; "/ninegurmukhi", [0x0A6F]; "/ninehackarabic", [0x0669];
"/ninehangzhou", [0x3029]; "/nineideographicparen", [0x3228]; "/nineinferior",
[0x2089]; "/ninemonospace", [0xFF19]; "/nineoldstyle", [0xF739]; "/nineparen",
[0x247C]; "/nineperiod", [0x2490]; "/ninepersian", [0x06F9]; "/nineroman",
[0x2178]; "/ninesuperior", [0x2079]; "/nineteencircle", [0x2472];
"/nineteenparen", [0x2486]; "/nineteenperiod", [0x249A]; "/ninethai", [0x0E59];
"/nj", [0x01CC]; "/njecyrillic", [0x045A]; "/nkatakana", [0x30F3];
"/nkatakanahalfwidth", [0xFF9D]; "/nlegrightlong", [0x019E]; "/nlinebelow",
[0x1E49]; "/nmonospace", [0xFF4E]; "/nmsquare", [0x339A]; "/nnabengali",
[0x09A3]; "/nnadeva", [0x0923]; "/nnagujarati", [0x0AA3]; "/nnagurmukhi",
[0x0A23]; "/nnnadeva", [0x0929]; "/nohiragana", [0x306E]; "/nokatakana",
[0x30CE]; "/nokatakanahalfwidth", [0xFF89]; "/nonbreakingspace", [0x00A0];
"/nonenthai", [0x0E13]; "/nonuthai", [0x0E19]; "/noonarabic", [0x0646];
"/noonfinalarabic", [0xFEE6]; "/noonghunnaarabic", [0x06BA];
"/noonghunnafinalarabic", [0xFB9F]; "/noonhehinitialarabic", [0xFEE7; 0xFEEC];
"/nooninitialarabic", [0xFEE7]; "/noonjeeminitialarabic", [0xFCD2];
"/noonjeemisolatedarabic", [0xFC4B]; "/noonmedialarabic", [0xFEE8];
"/noonmeeminitialarabic", [0xFCD5]; "/noonmeemisolatedarabic", [0xFC4E];
"/noonnoonfinalarabic", [0xFC8D]; "/notcontains", [0x220C]; "/notelement",
[0x2209]; "/notelementof", [0x2209]; "/notequal", [0x2260]; "/notgreater",
[0x226F]; "/notgreaternorequal", [0x2271]; "/notgreaternorless", [0x2279];
"/notidentical", [0x2262]; "/notless", [0x226E]; "/notlessnorequal", [0x2270];
"/notparallel", [0x2226]; "/notprecedes", [0x2280]; "/notsubset", [0x2284];
"/notsucceeds", [0x2281]; "/notsuperset", [0x2285]; "/nowarmenian", [0x0576];
"/nparen", [0x24A9]; "/nssquare", [0x33B1]; "/nsuperior", [0x207F]; "/ntilde",
[0x00F1]; "/nu", [0x03BD]; "/nuhiragana", [0x306C]; "/nukatakana", [0x30CC];
"/nukatakanahalfwidth", [0xFF87]; "/nuktabengali", [0x09BC]; "/nuktadeva",
[0x093C]; "/nuktagujarati", [0x0ABC]; "/nuktagurmukhi", [0x0A3C]; "/numbersign",
[0x0023]; "/numbersignmonospace", [0xFF03]; "/numbersignsmall", [0xFE5F];
"/numeralsigngreek", [0x0374]; "/numeralsignlowergreek", [0x0375]; "/numero",
[0x2116]; "/nun", [0x05E0]; "/nundagesh", [0xFB40]; "/nundageshhebrew",
[0xFB40]; "/nunhebrew", [0x05E0]; "/nvsquare", [0x33B5]; "/nwsquare", [0x33BB];
"/nyabengali", [0x099E]; "/nyadeva", [0x091E]; "/nyagujarati", [0x0A9E];
"/nyagurmukhi", [0x0A1E]; "/o", [0x006F]; "/oacute", [0x00F3]; "/oangthai",
[0x0E2D]; "/obarred", [0x0275]; "/obarredcyrillic", [0x04E9];
"/obarreddieresiscyrillic", [0x04EB]; "/obengali", [0x0993]; "/obopomofo",
[0x311B]; "/obreve", [0x014F]; "/ocandradeva", [0x0911]; "/ocandragujarati",
[0x0A91]; "/ocandravowelsigndeva", [0x0949]; "/ocandravowelsigngujarati",
[0x0AC9]; "/ocaron", [0x01D2]; "/ocircle", [0x24DE]; "/ocircumflex", [0x00F4];
"/ocircumflexacute", [0x1ED1]; "/ocircumflexdotbelow", [0x1ED9];
"/ocircumflexgrave", [0x1ED3]; "/ocircumflexhookabove", [0x1ED5];
"/ocircumflextilde", [0x1ED7]; "/ocyrillic", [0x043E]; "/odblacute", [0x0151];
"/odblgrave", [0x020D]; "/odeva", [0x0913]; "/odieresis", [0x00F6];
"/odieresiscyrillic", [0x04E7]; "/odotbelow", [0x1ECD]; "/oe", [0x0153];
"/oekorean", [0x315A]; "/ogonek", [0x02DB]; "/ogonekcmb", [0x0328]; "/ograve",
[0x00F2]; "/ogujarati", [0x0A93]; "/oharmenian", [0x0585]; "/ohiragana",
[0x304A]; "/ohookabove", [0x1ECF]; "/ohorn", [0x01A1]; "/ohornacute", [0x1EDB];
"/ohorndotbelow", [0x1EE3]; "/ohorngrave", [0x1EDD]; "/ohornhookabove",
[0x1EDF]; "/ohorntilde", [0x1EE1]; "/ohungarumlaut", [0x0151]; "/oi", [0x01A3];
"/oinvertedbreve", [0x020F]; "/okatakana", [0x30AA]; "/okatakanahalfwidth",
[0xFF75]; "/okorean", [0x3157]; "/olehebrew", [0x05AB]; "/omacron", [0x014D];
"/omacronacute", [0x1E53]; "/omacrongrave", [0x1E51]; "/omdeva", [0x0950];
"/omega", [0x03C9]; "/omega1", [0x03D6]; "/omegacyrillic", [0x0461];
"/omegalatinclosed", [0x0277]; "/omegaroundcyrillic", [0x047B];
"/omegatitlocyrillic", [0x047D]; "/omegatonos", [0x03CE]; "/omgujarati",
[0x0AD0]; "/omicron", [0x03BF]; "/omicrontonos", [0x03CC]; "/omonospace",
[0xFF4F]; "/one", [0x0031]; "/onearabic", [0x0661]; "/onebengali", [0x09E7];
"/onecircle", [0x2460]; "/onecircleinversesansserif", [0x278A]; "/onedeva",
[0x0967]; "/onedotenleader", [0x2024]; "/oneeighth", [0x215B]; "/onefitted",
[0xF6DC]; "/onegujarati", [0x0AE7]; "/onegurmukhi", [0x0A67]; "/onehackarabic",
[0x0661]; "/onehalf", [0x00BD]; "/onehangzhou", [0x3021];
"/oneideographicparen", [0x3220]; "/oneinferior", [0x2081]; "/onemonospace",
[0xFF11]; "/onenumeratorbengali", [0x09F4]; "/oneoldstyle", [0xF731];
"/oneparen", [0x2474]; "/oneperiod", [0x2488]; "/onepersian", [0x06F1];
"/onequarter", [0x00BC]; "/oneroman", [0x2170]; "/onesuperior", [0x00B9];
"/onethai", [0x0E51]; "/onethird", [0x2153]; "/oogonek", [0x01EB];
"/oogonekmacron", [0x01ED]; "/oogurmukhi", [0x0A13]; "/oomatragurmukhi",
[0x0A4B]; "/oopen", [0x0254]; "/oparen", [0x24AA]; "/openbullet", [0x25E6];
"/option", [0x2325]; "/ordfeminine", [0x00AA]; "/ordmasculine", [0x00BA];
"/orthogonal", [0x221F]; "/oshortdeva", [0x0912]; "/oshortvowelsigndeva",
[0x094A]; "/oslash", [0x00F8]; "/oslashacute", [0x01FF]; "/osmallhiragana",
[0x3049]; "/osmallkatakana", [0x30A9]; "/osmallkatakanahalfwidth", [0xFF6B];
"/ostrokeacute", [0x01FF]; "/osuperior", [0xF6F0]; "/otcyrillic", [0x047F];
"/otilde", [0x00F5]; "/otildeacute", [0x1E4D]; "/otildedieresis", [0x1E4F];
"/oubopomofo", [0x3121]; "/overline", [0x203E]; "/overlinecenterline", [0xFE4A];
"/overlinecmb", [0x0305]; "/overlinedashed", [0xFE49]; "/overlinedblwavy",
[0xFE4C]; "/overlinewavy", [0xFE4B]; "/overscore", [0x00AF];
"/ovowelsignbengali", [0x09CB]; "/ovowelsigndeva", [0x094B];
"/ovowelsigngujarati", [0x0ACB]; "/p", [0x0070]; "/paampssquare", [0x3380];
"/paasentosquare", [0x332B]; "/pabengali", [0x09AA]; "/pacute", [0x1E55];
"/padeva", [0x092A]; "/pagedown", [0x21DF]; "/pageup", [0x21DE]; "/pagujarati",
[0x0AAA]; "/pagurmukhi", [0x0A2A]; "/pahiragana", [0x3071]; "/paiyannoithai",
[0x0E2F]; "/pakatakana", [0x30D1]; "/palatalizationcyrilliccmb", [0x0484];
"/palochkacyrillic", [0x04C0]; "/pansioskorean", [0x317F]; "/paragraph",
[0x00B6]; "/parallel", [0x2225]; "/parenleft", [0x0028];
"/parenleftaltonearabic", [0xFD3E]; "/parenleftbt", [0xF8ED]; "/parenleftex",
[0xF8EC]; "/parenleftinferior", [0x208D]; "/parenleftmonospace", [0xFF08];
"/parenleftsmall", [0xFE59]; "/parenleftsuperior", [0x207D]; "/parenlefttp",
[0xF8EB]; "/parenleftvertical", [0xFE35]; "/parenright", [0x0029];
"/parenrightaltonearabic", [0xFD3F]; "/parenrightbt", [0xF8F8]; "/parenrightex",
[0xF8F7]; "/parenrightinferior", [0x208E]; "/parenrightmonospace", [0xFF09];
"/parenrightsmall", [0xFE5A]; "/parenrightsuperior", [0x207E]; "/parenrighttp",
[0xF8F6]; "/parenrightvertical", [0xFE36]; "/partialdiff", [0x2202];
"/paseqhebrew", [0x05C0]; "/pashtahebrew", [0x0599]; "/pasquare", [0x33A9];
"/patah", [0x05B7]; "/patah11", [0x05B7]; "/patah1d", [0x05B7]; "/patah2a",
[0x05B7]; "/patahhebrew", [0x05B7]; "/patahnarrowhebrew", [0x05B7];
"/patahquarterhebrew", [0x05B7]; "/patahwidehebrew", [0x05B7]; "/pazerhebrew",
[0x05A1]; "/pbopomofo", [0x3106]; "/pcircle", [0x24DF]; "/pdotaccent", [0x1E57];
"/pe", [0x05E4]; "/pecyrillic", [0x043F]; "/pedagesh", [0xFB44];
"/pedageshhebrew", [0xFB44]; "/peezisquare", [0x333B]; "/pefinaldageshhebrew",
[0xFB43]; "/peharabic", [0x067E]; "/peharmenian", [0x057A]; "/pehebrew",
[0x05E4]; "/pehfinalarabic", [0xFB57]; "/pehinitialarabic", [0xFB58];
"/pehiragana", [0x307A]; "/pehmedialarabic", [0xFB59]; "/pekatakana", [0x30DA];
"/pemiddlehookcyrillic", [0x04A7]; "/perafehebrew", [0xFB4E]; "/percent",
[0x0025]; "/percentarabic", [0x066A]; "/percentmonospace", [0xFF05];
"/percentsmall", [0xFE6A]; "/period", [0x002E]; "/periodarmenian", [0x0589];
"/periodcentered", [0x00B7]; "/periodhalfwidth", [0xFF61]; "/periodinferior",
[0xF6E7]; "/periodmonospace", [0xFF0E]; "/periodsmall", [0xFE52];
"/periodsuperior", [0xF6E8]; "/perispomenigreekcmb", [0x0342]; "/perpendicular",
[0x22A5]; "/perthousand", [0x2030]; "/peseta", [0x20A7]; "/pfsquare", [0x338A];
"/phabengali", [0x09AB]; "/phadeva", [0x092B]; "/phagujarati", [0x0AAB];
"/phagurmukhi", [0x0A2B]; "/phi", [0x03C6]; "/phi1", [0x03D5];
"/phieuphacirclekorean", [0x327A]; "/phieuphaparenkorean", [0x321A];
"/phieuphcirclekorean", [0x326C]; "/phieuphkorean", [0x314D];
"/phieuphparenkorean", [0x320C]; "/philatin", [0x0278]; "/phinthuthai",
[0x0E3A]; "/phisymbolgreek", [0x03D5]; "/phook", [0x01A5]; "/phophanthai",
[0x0E1E]; "/phophungthai", [0x0E1C]; "/phosamphaothai", [0x0E20]; "/pi",
[0x03C0]; "/pieupacirclekorean", [0x3273]; "/pieupaparenkorean", [0x3213];
"/pieupcieuckorean", [0x3176]; "/pieupcirclekorean", [0x3265];
"/pieupkiyeokkorean", [0x3172]; "/pieupkorean", [0x3142]; "/pieupparenkorean",
[0x3205]; "/pieupsioskiyeokkorean", [0x3174]; "/pieupsioskorean", [0x3144];
"/pieupsiostikeutkorean", [0x3175]; "/pieupthieuthkorean", [0x3177];
"/pieuptikeutkorean", [0x3173]; "/pihiragana", [0x3074]; "/pikatakana",
[0x30D4]; "/pisymbolgreek", [0x03D6]; "/piwrarmenian", [0x0583]; "/plus",
[0x002B]; "/plusbelowcmb", [0x031F]; "/pluscircle", [0x2295]; "/plusminus",
[0x00B1]; "/plusmod", [0x02D6]; "/plusmonospace", [0xFF0B]; "/plussmall",
[0xFE62]; "/plussuperior", [0x207A]; "/pmonospace", [0xFF50]; "/pmsquare",
[0x33D8]; "/pohiragana", [0x307D]; "/pointingindexdownwhite", [0x261F];
"/pointingindexleftwhite", [0x261C]; "/pointingindexrightwhite", [0x261E];
"/pointingindexupwhite", [0x261D]; "/pokatakana", [0x30DD]; "/poplathai",
[0x0E1B]; "/postalmark", [0x3012]; "/postalmarkface", [0x3020]; "/pparen",
[0x24AB]; "/precedes", [0x227A]; "/prescription", [0x211E]; "/primemod",
[0x02B9]; "/primereversed", [0x2035]; "/product", [0x220F]; "/projective",
[0x2305]; "/prolongedkana", [0x30FC]; "/propellor", [0x2318]; "/propersubset",
[0x2282]; "/propersuperset", [0x2283]; "/proportion", [0x2237]; "/proportional",
[0x221D]; "/psi", [0x03C8]; "/psicyrillic", [0x0471];
"/psilipneumatacyrilliccmb", [0x0486]; "/pssquare", [0x33B0]; "/puhiragana",
[0x3077]; "/pukatakana", [0x30D7]; "/pvsquare", [0x33B4]; "/pwsquare", [0x33BA];
"/q", [0x0071]; "/qadeva", [0x0958]; "/qadmahebrew", [0x05A8]; "/qafarabic",
[0x0642]; "/qaffinalarabic", [0xFED6]; "/qafinitialarabic", [0xFED7];
"/qafmedialarabic", [0xFED8]; "/qamats", [0x05B8]; "/qamats10", [0x05B8];
"/qamats1a", [0x05B8]; "/qamats1c", [0x05B8]; "/qamats27", [0x05B8];
"/qamats29", [0x05B8]; "/qamats33", [0x05B8]; "/qamatsde", [0x05B8];
"/qamatshebrew", [0x05B8]; "/qamatsnarrowhebrew", [0x05B8];
"/qamatsqatanhebrew", [0x05B8]; "/qamatsqatannarrowhebrew", [0x05B8];
"/qamatsqatanquarterhebrew", [0x05B8]; "/qamatsqatanwidehebrew", [0x05B8];
"/qamatsquarterhebrew", [0x05B8]; "/qamatswidehebrew", [0x05B8];
"/qarneyparahebrew", [0x059F]; "/qbopomofo", [0x3111]; "/qcircle", [0x24E0];
"/qhook", [0x02A0]; "/qmonospace", [0xFF51]; "/qof", [0x05E7]; "/qofdagesh",
[0xFB47]; "/qofdageshhebrew", [0xFB47]; "/qofhatafpatah", [0x05E7; 0x05B2];
"/qofhatafpatahhebrew", [0x05E7; 0x05B2]; "/qofhatafsegol", [0x05E7; 0x05B1];
"/qofhatafsegolhebrew", [0x05E7; 0x05B1]; "/qofhebrew", [0x05E7]; "/qofhiriq",
[0x05E7; 0x05B4]; "/qofhiriqhebrew", [0x05E7; 0x05B4]; "/qofholam", [0x05E7; 0x05B9];
"/qofholamhebrew", [0x05E7; 0x05B9]; "/qofpatah", [0x05E7; 0x05B7]; "/qofpatahhebrew",
[0x05E7; 0x05B7]; "/qofqamats", [0x05E7; 0x05B8]; "/qofqamatshebrew", [0x05E7;
0x05B8];
"/qofqubuts", [0x05E7; 0x05BB]; "/qofqubutshebrew", [0x05E7; 0x05BB]; "/qofsegol",
[0x05E7; 0x05B6]; "/qofsegolhebrew", [0x05E7; 0x05B6]; "/qofsheva", [0x05E7; 0x05B0];
"/qofshevahebrew", [0x05E7; 0x05B0]; "/qoftsere", [0x05E7; 0x05B5]; "/qoftserehebrew",
[0x05E7; 0x05B5]; "/qparen", [0x24AC]; "/quarternote", [0x2669]; "/qubuts",
[0x05BB]; "/qubuts18", [0x05BB]; "/qubuts25", [0x05BB]; "/qubuts31", [0x05BB];
"/qubutshebrew", [0x05BB]; "/qubutsnarrowhebrew", [0x05BB];
"/qubutsquarterhebrew", [0x05BB]; "/qubutswidehebrew", [0x05BB]; "/question",
[0x003F]; "/questionarabic", [0x061F]; "/questionarmenian", [0x055E];
"/questiondown", [0x00BF]; "/questiondownsmall", [0xF7BF]; "/questiongreek",
[0x037E]; "/questionmonospace", [0xFF1F]; "/questionsmall", [0xF73F];
"/quotedbl", [0x0022]; "/quotedblbase", [0x201E]; "/quotedblleft", [0x201C];
"/quotedblmonospace", [0xFF02]; "/quotedblprime", [0x301E];
"/quotedblprimereversed", [0x301D]; "/quotedblright", [0x201D]; "/quoteleft",
[0x2018]; "/quoteleftreversed", [0x201B]; "/quotereversed", [0x201B];
"/quoteright", [0x2019]; "/quoterightn", [0x0149]; "/quotesinglbase", [0x201A];
"/quotesingle", [0x0027]; "/quotesinglemonospace", [0xFF07]; "/r", [0x0072];
"/raarmenian", [0x057C]; "/rabengali", [0x09B0]; "/racute", [0x0155]; "/radeva",
[0x0930]; "/radical", [0x221A]; "/radicalex", [0xF8E5]; "/radoverssquare",
[0x33AE]; "/radoverssquaredsquare", [0x33AF]; "/radsquare", [0x33AD]; "/rafe",
[0x05BF]; "/rafehebrew", [0x05BF]; "/ragujarati", [0x0AB0]; "/ragurmukhi",
[0x0A30]; "/rahiragana", [0x3089]; "/rakatakana", [0x30E9];
"/rakatakanahalfwidth", [0xFF97]; "/ralowerdiagonalbengali", [0x09F1];
"/ramiddlediagonalbengali", [0x09F0]; "/ramshorn", [0x0264]; "/ratio", [0x2236];
"/rbopomofo", [0x3116]; "/rcaron", [0x0159]; "/rcedilla", [0x0157]; "/rcircle",
[0x24E1]; "/rcommaaccent", [0x0157]; "/rdblgrave", [0x0211]; "/rdotaccent",
[0x1E59]; "/rdotbelow", [0x1E5B]; "/rdotbelowmacron", [0x1E5D];
"/referencemark", [0x203B]; "/reflexsubset", [0x2286]; "/reflexsuperset",
[0x2287]; "/registered", [0x00AE]; "/registersans", [0xF8E8]; "/registerserif",
[0xF6DA]; "/reharabic", [0x0631]; "/reharmenian", [0x0580]; "/rehfinalarabic",
[0xFEAE]; "/rehiragana", [0x308C]; "/rehyehaleflamarabic", [0x0631; 0xFEF3;
0xFE8E;
 0x0644]; "/rekatakana", [0x30EC]; "/rekatakanahalfwidth", [0xFF9A]; "/resh",
[0x05E8]; "/reshdageshhebrew", [0xFB48]; "/reshhatafpatah", [0x05E8; 0x05B2];
"/reshhatafpatahhebrew", [0x05E8; 0x05B2]; "/reshhatafsegol", [0x05E8; 0x05B1];
"/reshhatafsegolhebrew", [0x05E8; 0x05B1]; "/reshhebrew", [0x05E8]; "/reshhiriq",
[0x05E8; 0x05B4]; "/reshhiriqhebrew", [0x05E8; 0x05B4]; "/reshholam", [0x05E8; 0x05B9];
"/reshholamhebrew", [0x05E8; 0x05B9]; "/reshpatah", [0x05E8; 0x05B7];
"/reshpatahhebrew", [0x05E8; 0x05B7]; "/reshqamats", [0x05E8; 0x05B8];
"/reshqamatshebrew", [0x05E8; 0x05B8]; "/reshqubuts", [0x05E8; 0x05BB];
"/reshqubutshebrew", [0x05E8; 0x05BB]; "/reshsegol", [0x05E8; 0x05B6];
"/reshsegolhebrew", [0x05E8; 0x05B6]; "/reshsheva", [0x05E8; 0x05B0];
"/reshshevahebrew", [0x05E8; 0x05B0]; "/reshtsere", [0x05E8; 0x05B5];
"/reshtserehebrew", [0x05E8; 0x05B5]; "/reversedtilde", [0x223D]; "/reviahebrew",
[0x0597]; "/reviamugrashhebrew", [0x0597]; "/revlogicalnot", [0x2310];
"/rfishhook", [0x027E]; "/rfishhookreversed", [0x027F]; "/rhabengali", [0x09DD];
"/rhadeva", [0x095D]; "/rho", [0x03C1]; "/rhook", [0x027D]; "/rhookturned",
[0x027B]; "/rhookturnedsuperior", [0x02B5]; "/rhosymbolgreek", [0x03F1];
"/rhotichookmod", [0x02DE]; "/rieulacirclekorean", [0x3271];
"/rieulaparenkorean", [0x3211]; "/rieulcirclekorean", [0x3263];
"/rieulhieuhkorean", [0x3140]; "/rieulkiyeokkorean", [0x313A];
"/rieulkiyeoksioskorean", [0x3169]; "/rieulkorean", [0x3139];
"/rieulmieumkorean", [0x313B]; "/rieulpansioskorean", [0x316C];
"/rieulparenkorean", [0x3203]; "/rieulphieuphkorean", [0x313F];
"/rieulpieupkorean", [0x313C]; "/rieulpieupsioskorean", [0x316B];
"/rieulsioskorean", [0x313D]; "/rieulthieuthkorean", [0x313E];
"/rieultikeutkorean", [0x316A]; "/rieulyeorinhieuhkorean", [0x316D];
"/rightangle", [0x221F]; "/righttackbelowcmb", [0x0319]; "/righttriangle",
[0x22BF]; "/rihiragana", [0x308A]; "/rikatakana", [0x30EA];
"/rikatakanahalfwidth", [0xFF98]; "/ring", [0x02DA]; "/ringbelowcmb", [0x0325];
"/ringcmb", [0x030A]; "/ringhalfleft", [0x02BF]; "/ringhalfleftarmenian",
[0x0559]; "/ringhalfleftbelowcmb", [0x031C]; "/ringhalfleftcentered", [0x02D3];
"/ringhalfright", [0x02BE]; "/ringhalfrightbelowcmb", [0x0339];
"/ringhalfrightcentered", [0x02D2]; "/rinvertedbreve", [0x0213];
"/rittorusquare", [0x3351]; "/rlinebelow", [0x1E5F]; "/rlongleg", [0x027C];
"/rlonglegturned", [0x027A]; "/rmonospace", [0xFF52]; "/rohiragana", [0x308D];
"/rokatakana", [0x30ED]; "/rokatakanahalfwidth", [0xFF9B]; "/roruathai",
[0x0E23]; "/rparen", [0x24AD]; "/rrabengali", [0x09DC]; "/rradeva", [0x0931];
"/rragurmukhi", [0x0A5C]; "/rreharabic", [0x0691]; "/rrehfinalarabic", [0xFB8D];
"/rrvocalicbengali", [0x09E0]; "/rrvocalicdeva", [0x0960]; "/rrvocalicgujarati",
[0x0AE0]; "/rrvocalicvowelsignbengali", [0x09C4]; "/rrvocalicvowelsigndeva",
[0x0944]; "/rrvocalicvowelsigngujarati", [0x0AC4]; "/rsuperior", [0xF6F1];
"/rtblock", [0x2590]; "/rturned", [0x0279]; "/rturnedsuperior", [0x02B4];
"/ruhiragana", [0x308B]; "/rukatakana", [0x30EB]; "/rukatakanahalfwidth",
[0xFF99]; "/rupeemarkbengali", [0x09F2]; "/rupeesignbengali", [0x09F3];
"/rupiah", [0xF6DD]; "/ruthai", [0x0E24]; "/rvocalicbengali", [0x098B];
"/rvocalicdeva", [0x090B]; "/rvocalicgujarati", [0x0A8B];
"/rvocalicvowelsignbengali", [0x09C3]; "/rvocalicvowelsigndeva", [0x0943];
"/rvocalicvowelsigngujarati", [0x0AC3]; "/s", [0x0073]; "/sabengali", [0x09B8];
"/sacute", [0x015B]; "/sacutedotaccent", [0x1E65]; "/sadarabic", [0x0635];
"/sadeva", [0x0938]; "/sadfinalarabic", [0xFEBA]; "/sadinitialarabic", [0xFEBB];
"/sadmedialarabic", [0xFEBC]; "/sagujarati", [0x0AB8]; "/sagurmukhi", [0x0A38];
"/sahiragana", [0x3055]; "/sakatakana", [0x30B5]; "/sakatakanahalfwidth",
[0xFF7B]; "/sallallahoualayhewasallamarabic", [0xFDFA]; "/samekh", [0x05E1];
"/samekhdagesh", [0xFB41]; "/samekhdageshhebrew", [0xFB41]; "/samekhhebrew",
[0x05E1]; "/saraaathai", [0x0E32]; "/saraaethai", [0x0E41];
"/saraaimaimalaithai", [0x0E44]; "/saraaimaimuanthai", [0x0E43]; "/saraamthai",
[0x0E33]; "/saraathai", [0x0E30]; "/saraethai", [0x0E40]; "/saraiileftthai",
[0xF886]; "/saraiithai", [0x0E35]; "/saraileftthai", [0xF885]; "/saraithai",
[0x0E34]; "/saraothai", [0x0E42]; "/saraueeleftthai", [0xF888]; "/saraueethai",
[0x0E37]; "/saraueleftthai", [0xF887]; "/sarauethai", [0x0E36]; "/sarauthai",
[0x0E38]; "/sarauuthai", [0x0E39]; "/sbopomofo", [0x3119]; "/scaron", [0x0161];
"/scarondotaccent", [0x1E67]; "/scedilla", [0x015F]; "/schwa", [0x0259];
"/schwacyrillic", [0x04D9]; "/schwadieresiscyrillic", [0x04DB]; "/schwahook",
[0x025A]; "/scircle", [0x24E2]; "/scircumflex", [0x015D]; "/scommaaccent",
[0x0219]; "/sdotaccent", [0x1E61]; "/sdotbelow", [0x1E63];
"/sdotbelowdotaccent", [0x1E69]; "/seagullbelowcmb", [0x033C]; "/second",
[0x2033]; "/secondtonechinese", [0x02CA]; "/section", [0x00A7]; "/seenarabic",
[0x0633]; "/seenfinalarabic", [0xFEB2]; "/seeninitialarabic", [0xFEB3];
"/seenmedialarabic", [0xFEB4]; "/segol", [0x05B6]; "/segol13", [0x05B6];
"/segol1f", [0x05B6]; "/segol2c", [0x05B6]; "/segolhebrew", [0x05B6];
"/segolnarrowhebrew", [0x05B6]; "/segolquarterhebrew", [0x05B6];
"/segoltahebrew", [0x0592]; "/segolwidehebrew", [0x05B6]; "/seharmenian",
[0x057D]; "/sehiragana", [0x305B]; "/sekatakana", [0x30BB];
"/sekatakanahalfwidth", [0xFF7E]; "/semicolon", [0x003B]; "/semicolonarabic",
[0x061B]; "/semicolonmonospace", [0xFF1B]; "/semicolonsmall", [0xFE54];
"/semivoicedmarkkana", [0x309C]; "/semivoicedmarkkanahalfwidth", [0xFF9F];
"/sentisquare", [0x3322]; "/sentosquare", [0x3323]; "/seven", [0x0037];
"/sevenarabic", [0x0667]; "/sevenbengali", [0x09ED]; "/sevencircle", [0x2466];
"/sevencircleinversesansserif", [0x2790]; "/sevendeva", [0x096D];
"/seveneighths", [0x215E]; "/sevengujarati", [0x0AED]; "/sevengurmukhi",
[0x0A6D]; "/sevenhackarabic", [0x0667]; "/sevenhangzhou", [0x3027];
"/sevenideographicparen", [0x3226]; "/seveninferior", [0x2087];
"/sevenmonospace", [0xFF17]; "/sevenoldstyle", [0xF737]; "/sevenparen",
[0x247A]; "/sevenperiod", [0x248E]; "/sevenpersian", [0x06F7]; "/sevenroman",
[0x2176]; "/sevensuperior", [0x2077]; "/seventeencircle", [0x2470];
"/seventeenparen", [0x2484]; "/seventeenperiod", [0x2498]; "/seventhai",
[0x0E57]; "/sfthyphen", [0x00AD]; "/shaarmenian", [0x0577]; "/shabengali",
[0x09B6]; "/shacyrillic", [0x0448]; "/shaddaarabic", [0x0651];
"/shaddadammaarabic", [0xFC61]; "/shaddadammatanarabic", [0xFC5E];
"/shaddafathaarabic", [0xFC60]; "/shaddafathatanarabic", [0x0651; 0x064B];
"/shaddakasraarabic", [0xFC62]; "/shaddakasratanarabic", [0xFC5F]; "/shade",
[0x2592]; "/shadedark", [0x2593]; "/shadelight", [0x2591]; "/shademedium",
[0x2592]; "/shadeva", [0x0936]; "/shagujarati", [0x0AB6]; "/shagurmukhi",
[0x0A36]; "/shalshelethebrew", [0x0593]; "/shbopomofo", [0x3115];
"/shchacyrillic", [0x0449]; "/sheenarabic", [0x0634]; "/sheenfinalarabic",
[0xFEB6]; "/sheeninitialarabic", [0xFEB7]; "/sheenmedialarabic", [0xFEB8];
"/sheicoptic", [0x03E3]; "/sheqel", [0x20AA]; "/sheqelhebrew", [0x20AA];
"/sheva", [0x05B0]; "/sheva115", [0x05B0]; "/sheva15", [0x05B0]; "/sheva22",
[0x05B0]; "/sheva2e", [0x05B0]; "/shevahebrew", [0x05B0]; "/shevanarrowhebrew",
[0x05B0]; "/shevaquarterhebrew", [0x05B0]; "/shevawidehebrew", [0x05B0];
"/shhacyrillic", [0x04BB]; "/shimacoptic", [0x03ED]; "/shin", [0x05E9];
"/shindagesh", [0xFB49]; "/shindageshhebrew", [0xFB49]; "/shindageshshindot",
[0xFB2C]; "/shindageshshindothebrew", [0xFB2C]; "/shindageshsindot", [0xFB2D];
"/shindageshsindothebrew", [0xFB2D]; "/shindothebrew", [0x05C1]; "/shinhebrew",
[0x05E9]; "/shinshindot", [0xFB2A]; "/shinshindothebrew", [0xFB2A];
"/shinsindot", [0xFB2B]; "/shinsindothebrew", [0xFB2B]; "/shook", [0x0282];
"/sigma", [0x03C3]; "/sigma1", [0x03C2]; "/sigmafinal", [0x03C2];
"/sigmalunatesymbolgreek", [0x03F2]; "/sihiragana", [0x3057]; "/sikatakana",
[0x30B7]; "/sikatakanahalfwidth", [0xFF7C]; "/siluqhebrew", [0x05BD];
"/siluqlefthebrew", [0x05BD]; "/similar", [0x223C]; "/sindothebrew", [0x05C2];
"/siosacirclekorean", [0x3274]; "/siosaparenkorean", [0x3214];
"/sioscieuckorean", [0x317E]; "/sioscirclekorean", [0x3266];
"/sioskiyeokkorean", [0x317A]; "/sioskorean", [0x3145]; "/siosnieunkorean",
[0x317B]; "/siosparenkorean", [0x3206]; "/siospieupkorean", [0x317D];
"/siostikeutkorean", [0x317C]; "/six", [0x0036]; "/sixarabic", [0x0666];
"/sixbengali", [0x09EC]; "/sixcircle", [0x2465]; "/sixcircleinversesansserif",
[0x278F]; "/sixdeva", [0x096C]; "/sixgujarati", [0x0AEC]; "/sixgurmukhi",
[0x0A6C]; "/sixhackarabic", [0x0666]; "/sixhangzhou", [0x3026];
"/sixideographicparen", [0x3225]; "/sixinferior", [0x2086]; "/sixmonospace",
[0xFF16]; "/sixoldstyle", [0xF736]; "/sixparen", [0x2479]; "/sixperiod",
[0x248D]; "/sixpersian", [0x06F6]; "/sixroman", [0x2175]; "/sixsuperior",
[0x2076]; "/sixteencircle", [0x246F]; "/sixteencurrencydenominatorbengali",
[0x09F9]; "/sixteenparen", [0x2483]; "/sixteenperiod", [0x2497]; "/sixthai",
[0x0E56]; "/slash", [0x002F]; "/slashmonospace", [0xFF0F]; "/slong", [0x017F];
"/slongdotaccent", [0x1E9B]; "/smileface", [0x263A]; "/smonospace", [0xFF53];
"/sofpasuqhebrew", [0x05C3]; "/softhyphen", [0x00AD]; "/softsigncyrillic",
[0x044C]; "/sohiragana", [0x305D]; "/sokatakana", [0x30BD];
"/sokatakanahalfwidth", [0xFF7F]; "/soliduslongoverlaycmb", [0x0338];
"/solidusshortoverlaycmb", [0x0337]; "/sorusithai", [0x0E29]; "/sosalathai",
[0x0E28]; "/sosothai", [0x0E0B]; "/sosuathai", [0x0E2A]; "/space", [0x0020];
"/spacehackarabic", [0x0020]; "/spade", [0x2660]; "/spadesuitblack", [0x2660];
"/spadesuitwhite", [0x2664]; "/sparen", [0x24AE]; "/squarebelowcmb", [0x033B];
"/squarecc", [0x33C4]; "/squarecm", [0x339D]; "/squarediagonalcrosshatchfill",
[0x25A9]; "/squarehorizontalfill", [0x25A4]; "/squarekg", [0x338F]; "/squarekm",
[0x339E]; "/squarekmcapital", [0x33CE]; "/squareln", [0x33D1]; "/squarelog",
[0x33D2]; "/squaremg", [0x338E]; "/squaremil", [0x33D5]; "/squaremm", [0x339C];
"/squaremsquared", [0x33A1]; "/squareorthogonalcrosshatchfill", [0x25A6];
"/squareupperlefttolowerrightfill", [0x25A7];
"/squareupperrighttolowerleftfill", [0x25A8]; "/squareverticalfill", [0x25A5];
"/squarewhitewithsmallblack", [0x25A3]; "/srsquare", [0x33DB]; "/ssabengali",
[0x09B7]; "/ssadeva", [0x0937]; "/ssagujarati", [0x0AB7]; "/ssangcieuckorean",
[0x3149]; "/ssanghieuhkorean", [0x3185]; "/ssangieungkorean", [0x3180];
"/ssangkiyeokkorean", [0x3132]; "/ssangnieunkorean", [0x3165];
"/ssangpieupkorean", [0x3143]; "/ssangsioskorean", [0x3146];
"/ssangtikeutkorean", [0x3138]; "/ssuperior", [0xF6F2]; "/sterling", [0x00A3];
"/sterlingmonospace", [0xFFE1]; "/strokelongoverlaycmb", [0x0336];
"/strokeshortoverlaycmb", [0x0335]; "/subset", [0x2282]; "/subsetnotequal",
[0x228A]; "/subsetorequal", [0x2286]; "/succeeds", [0x227B]; "/suchthat",
[0x220B]; "/suhiragana", [0x3059]; "/sukatakana", [0x30B9];
"/sukatakanahalfwidth", [0xFF7D]; "/sukunarabic", [0x0652]; "/summation",
[0x2211]; "/sun", [0x263C]; "/superset", [0x2283]; "/supersetnotequal",
[0x228B]; "/supersetorequal", [0x2287]; "/svsquare", [0x33DC];
"/syouwaerasquare", [0x337C]; "/t", [0x0074]; "/tabengali", [0x09A4];
"/tackdown", [0x22A4]; "/tackleft", [0x22A3]; "/tadeva", [0x0924];
"/tagujarati", [0x0AA4]; "/tagurmukhi", [0x0A24]; "/taharabic", [0x0637];
"/tahfinalarabic", [0xFEC2]; "/tahinitialarabic", [0xFEC3]; "/tahiragana",
[0x305F]; "/tahmedialarabic", [0xFEC4]; "/taisyouerasquare", [0x337D];
"/takatakana", [0x30BF]; "/takatakanahalfwidth", [0xFF80]; "/tatweelarabic",
[0x0640]; "/tau", [0x03C4]; "/tav", [0x05EA]; "/tavdages", [0xFB4A];
"/tavdagesh", [0xFB4A]; "/tavdageshhebrew", [0xFB4A]; "/tavhebrew", [0x05EA];
"/tbar", [0x0167]; "/tbopomofo", [0x310A]; "/tcaron", [0x0165]; "/tccurl",
[0x02A8]; "/tcedilla", [0x0163]; "/tcheharabic", [0x0686]; "/tchehfinalarabic",
[0xFB7B]; "/tchehinitialarabic", [0xFB7C]; "/tchehmedialarabic", [0xFB7D];
"/tchehmeeminitialarabic", [0xFB7C; 0xFEE4]; "/tcircle", [0x24E3];
"/tcircumflexbelow", [0x1E71]; "/tcommaaccent", [0x0163]; "/tdieresis",
[0x1E97]; "/tdotaccent", [0x1E6B]; "/tdotbelow", [0x1E6D]; "/tecyrillic",
[0x0442]; "/tedescendercyrillic", [0x04AD]; "/teharabic", [0x062A];
"/tehfinalarabic", [0xFE96]; "/tehhahinitialarabic", [0xFCA2];
"/tehhahisolatedarabic", [0xFC0C]; "/tehinitialarabic", [0xFE97]; "/tehiragana",
[0x3066]; "/tehjeeminitialarabic", [0xFCA1]; "/tehjeemisolatedarabic", [0xFC0B];
"/tehmarbutaarabic", [0x0629]; "/tehmarbutafinalarabic", [0xFE94];
"/tehmedialarabic", [0xFE98]; "/tehmeeminitialarabic", [0xFCA4];
"/tehmeemisolatedarabic", [0xFC0E]; "/tehnoonfinalarabic", [0xFC73];
"/tekatakana", [0x30C6]; "/tekatakanahalfwidth", [0xFF83]; "/telephone",
[0x2121]; "/telephoneblack", [0x260E]; "/telishagedolahebrew", [0x05A0];
"/telishaqetanahebrew", [0x05A9]; "/tencircle", [0x2469];
"/tenideographicparen", [0x3229]; "/tenparen", [0x247D]; "/tenperiod", [0x2491];
"/tenroman", [0x2179]; "/tesh", [0x02A7]; "/tet", [0x05D8]; "/tetdagesh",
[0xFB38]; "/tetdageshhebrew", [0xFB38]; "/tethebrew", [0x05D8];
"/tetsecyrillic", [0x04B5]; "/tevirhebrew", [0x059B]; "/tevirlefthebrew",
[0x059B]; "/thabengali", [0x09A5]; "/thadeva", [0x0925]; "/thagujarati",
[0x0AA5]; "/thagurmukhi", [0x0A25]; "/thalarabic", [0x0630]; "/thalfinalarabic",
[0xFEAC]; "/thanthakhatlowleftthai", [0xF898]; "/thanthakhatlowrightthai",
[0xF897]; "/thanthakhatthai", [0x0E4C]; "/thanthakhatupperleftthai", [0xF896];
"/theharabic", [0x062B]; "/thehfinalarabic", [0xFE9A]; "/thehinitialarabic",
[0xFE9B]; "/thehmedialarabic", [0xFE9C]; "/thereexists", [0x2203]; "/therefore",
[0x2234]; "/theta", [0x03B8]; "/theta1", [0x03D1]; "/thetasymbolgreek",
[0x03D1]; "/thieuthacirclekorean", [0x3279]; "/thieuthaparenkorean", [0x3219];
"/thieuthcirclekorean", [0x326B]; "/thieuthkorean", [0x314C];
"/thieuthparenkorean", [0x320B]; "/thirteencircle", [0x246C]; "/thirteenparen",
[0x2480]; "/thirteenperiod", [0x2494]; "/thonangmonthothai", [0x0E11]; "/thook",
[0x01AD]; "/thophuthaothai", [0x0E12]; "/thorn", [0x00FE]; "/thothahanthai",
[0x0E17]; "/thothanthai", [0x0E10]; "/thothongthai", [0x0E18]; "/thothungthai",
[0x0E16]; "/thousandcyrillic", [0x0482]; "/thousandsseparatorarabic", [0x066C];
"/thousandsseparatorpersian", [0x066C]; "/three", [0x0033]; "/threearabic",
[0x0663]; "/threebengali", [0x09E9]; "/threecircle", [0x2462];
"/threecircleinversesansserif", [0x278C]; "/threedeva", [0x0969];
"/threeeighths", [0x215C]; "/threegujarati", [0x0AE9]; "/threegurmukhi",
[0x0A69]; "/threehackarabic", [0x0663]; "/threehangzhou", [0x3023];
"/threeideographicparen", [0x3222]; "/threeinferior", [0x2083];
"/threemonospace", [0xFF13]; "/threenumeratorbengali", [0x09F6];
"/threeoldstyle", [0xF733]; "/threeparen", [0x2476]; "/threeperiod", [0x248A];
"/threepersian", [0x06F3]; "/threequarters", [0x00BE]; "/threequartersemdash",
[0xF6DE]; "/threeroman", [0x2172]; "/threesuperior", [0x00B3]; "/threethai",
[0x0E53]; "/thzsquare", [0x3394]; "/tihiragana", [0x3061]; "/tikatakana",
[0x30C1]; "/tikatakanahalfwidth", [0xFF81]; "/tikeutacirclekorean", [0x3270];
"/tikeutaparenkorean", [0x3210]; "/tikeutcirclekorean", [0x3262];
"/tikeutkorean", [0x3137]; "/tikeutparenkorean", [0x3202]; "/tilde", [0x02DC];
"/tildebelowcmb", [0x0330]; "/tildecmb", [0x0303]; "/tildecomb", [0x0303];
"/tildedoublecmb", [0x0360]; "/tildeoperator", [0x223C]; "/tildeoverlaycmb",
[0x0334]; "/tildeverticalcmb", [0x033E]; "/timescircle", [0x2297];
"/tipehahebrew", [0x0596]; "/tipehalefthebrew", [0x0596]; "/tippigurmukhi",
[0x0A70]; "/titlocyrilliccmb", [0x0483]; "/tiwnarmenian", [0x057F];
"/tlinebelow", [0x1E6F]; "/tmonospace", [0xFF54]; "/toarmenian", [0x0569];
"/tohiragana", [0x3068]; "/tokatakana", [0x30C8]; "/tokatakanahalfwidth",
[0xFF84]; "/tonebarextrahighmod", [0x02E5]; "/tonebarextralowmod", [0x02E9];
"/tonebarhighmod", [0x02E6]; "/tonebarlowmod", [0x02E8]; "/tonebarmidmod",
[0x02E7]; "/tonefive", [0x01BD]; "/tonesix", [0x0185]; "/tonetwo", [0x01A8];
"/tonos", [0x0384]; "/tonsquare", [0x3327]; "/topatakthai", [0x0E0F];
"/tortoiseshellbracketleft", [0x3014]; "/tortoiseshellbracketleftsmall",
[0xFE5D]; "/tortoiseshellbracketleftvertical", [0xFE39];
"/tortoiseshellbracketright", [0x3015]; "/tortoiseshellbracketrightsmall",
[0xFE5E]; "/tortoiseshellbracketrightvertical", [0xFE3A]; "/totaothai",
[0x0E15]; "/tpalatalhook", [0x01AB]; "/tparen", [0x24AF]; "/trademark",
[0x2122]; "/trademarksans", [0xF8EA]; "/trademarkserif", [0xF6DB];
"/tretroflexhook", [0x0288]; "/triagdn", [0x25BC]; "/triaglf", [0x25C4];
"/triagrt", [0x25BA]; "/triagup", [0x25B2]; "/ts", [0x02A6]; "/tsadi", [0x05E6];
"/tsadidagesh", [0xFB46]; "/tsadidageshhebrew", [0xFB46]; "/tsadihebrew",
[0x05E6]; "/tsecyrillic", [0x0446]; "/tsere", [0x05B5]; "/tsere12", [0x05B5];
"/tsere1e", [0x05B5]; "/tsere2b", [0x05B5]; "/tserehebrew", [0x05B5];
"/tserenarrowhebrew", [0x05B5]; "/tserequarterhebrew", [0x05B5];
"/tserewidehebrew", [0x05B5]; "/tshecyrillic", [0x045B]; "/tsuperior", [0xF6F3];
"/ttabengali", [0x099F]; "/ttadeva", [0x091F]; "/ttagujarati", [0x0A9F];
"/ttagurmukhi", [0x0A1F]; "/tteharabic", [0x0679]; "/ttehfinalarabic", [0xFB67];
"/ttehinitialarabic", [0xFB68]; "/ttehmedialarabic", [0xFB69]; "/tthabengali",
[0x09A0]; "/tthadeva", [0x0920]; "/tthagujarati", [0x0AA0]; "/tthagurmukhi",
[0x0A20]; "/tturned", [0x0287]; "/tuhiragana", [0x3064]; "/tukatakana",
[0x30C4]; "/tukatakanahalfwidth", [0xFF82]; "/tusmallhiragana", [0x3063];
"/tusmallkatakana", [0x30C3]; "/tusmallkatakanahalfwidth", [0xFF6F];
"/twelvecircle", [0x246B]; "/twelveparen", [0x247F]; "/twelveperiod", [0x2493];
"/twelveroman", [0x217B]; "/twentycircle", [0x2473]; "/twentyhangzhou",
[0x5344]; "/twentyparen", [0x2487]; "/twentyperiod", [0x249B]; "/two", [0x0032];
"/twoarabic", [0x0662]; "/twobengali", [0x09E8]; "/twocircle", [0x2461];
"/twocircleinversesansserif", [0x278B]; "/twodeva", [0x0968]; "/twodotenleader",
[0x2025]; "/twodotleader", [0x2025]; "/twodotleadervertical", [0xFE30];
"/twogujarati", [0x0AE8]; "/twogurmukhi", [0x0A68]; "/twohackarabic", [0x0662];
"/twohangzhou", [0x3022]; "/twoideographicparen", [0x3221]; "/twoinferior",
[0x2082]; "/twomonospace", [0xFF12]; "/twonumeratorbengali", [0x09F5];
"/twooldstyle", [0xF732]; "/twoparen", [0x2475]; "/twoperiod", [0x2489];
"/twopersian", [0x06F2]; "/tworoman", [0x2171]; "/twostroke", [0x01BB];
"/twosuperior", [0x00B2]; "/twothai", [0x0E52]; "/twothirds", [0x2154]; "/u",
[0x0075]; "/uacute", [0x00FA]; "/ubar", [0x0289]; "/ubengali", [0x0989];
"/ubopomofo", [0x3128]; "/ubreve", [0x016D]; "/ucaron", [0x01D4]; "/ucircle",
[0x24E4]; "/ucircumflex", [0x00FB]; "/ucircumflexbelow", [0x1E77]; "/ucyrillic",
[0x0443]; "/udattadeva", [0x0951]; "/udblacute", [0x0171]; "/udblgrave",
[0x0215]; "/udeva", [0x0909]; "/udieresis", [0x00FC]; "/udieresisacute",
[0x01D8]; "/udieresisbelow", [0x1E73]; "/udieresiscaron", [0x01DA];
"/udieresiscyrillic", [0x04F1]; "/udieresisgrave", [0x01DC]; "/udieresismacron",
[0x01D6]; "/udotbelow", [0x1EE5]; "/ugrave", [0x00F9]; "/ugujarati", [0x0A89];
"/ugurmukhi", [0x0A09]; "/uhiragana", [0x3046]; "/uhookabove", [0x1EE7];
"/uhorn", [0x01B0]; "/uhornacute", [0x1EE9]; "/uhorndotbelow", [0x1EF1];
"/uhorngrave", [0x1EEB]; "/uhornhookabove", [0x1EED]; "/uhorntilde", [0x1EEF];
"/uhungarumlaut", [0x0171]; "/uhungarumlautcyrillic", [0x04F3];
"/uinvertedbreve", [0x0217]; "/ukatakana", [0x30A6]; "/ukatakanahalfwidth",
[0xFF73]; "/ukcyrillic", [0x0479]; "/ukorean", [0x315C]; "/umacron", [0x016B];
"/umacroncyrillic", [0x04EF]; "/umacrondieresis", [0x1E7B]; "/umatragurmukhi",
[0x0A41]; "/umonospace", [0xFF55]; "/underscore", [0x005F]; "/underscoredbl",
[0x2017]; "/underscoremonospace", [0xFF3F]; "/underscorevertical", [0xFE33];
"/underscorewavy", [0xFE4F]; "/union", [0x222A]; "/universal", [0x2200];
"/uogonek", [0x0173]; "/uparen", [0x24B0]; "/upblock", [0x2580];
"/upperdothebrew", [0x05C4]; "/upsilon", [0x03C5]; "/upsilondieresis", [0x03CB];
"/upsilondieresistonos", [0x03B0]; "/upsilonlatin", [0x028A]; "/upsilontonos",
[0x03CD]; "/uptackbelowcmb", [0x031D]; "/uptackmod", [0x02D4]; "/uragurmukhi",
[0x0A73]; "/uring", [0x016F]; "/ushortcyrillic", [0x045E]; "/usmallhiragana",
[0x3045]; "/usmallkatakana", [0x30A5]; "/usmallkatakanahalfwidth", [0xFF69];
"/ustraightcyrillic", [0x04AF]; "/ustraightstrokecyrillic", [0x04B1]; "/utilde",
[0x0169]; "/utildeacute", [0x1E79]; "/utildebelow", [0x1E75]; "/uubengali",
[0x098A]; "/uudeva", [0x090A]; "/uugujarati", [0x0A8A]; "/uugurmukhi", [0x0A0A];
"/uumatragurmukhi", [0x0A42]; "/uuvowelsignbengali", [0x09C2];
"/uuvowelsigndeva", [0x0942]; "/uuvowelsigngujarati", [0x0AC2];
"/uvowelsignbengali", [0x09C1]; "/uvowelsigndeva", [0x0941];
"/uvowelsigngujarati", [0x0AC1]; "/v", [0x0076]; "/vadeva", [0x0935];
"/vagujarati", [0x0AB5]; "/vagurmukhi", [0x0A35]; "/vakatakana", [0x30F7];
"/vav", [0x05D5]; "/vavdagesh", [0xFB35]; "/vavdagesh65", [0xFB35];
"/vavdageshhebrew", [0xFB35]; "/vavhebrew", [0x05D5]; "/vavholam", [0xFB4B];
"/vavholamhebrew", [0xFB4B]; "/vavvavhebrew", [0x05F0]; "/vavyodhebrew",
[0x05F1]; "/vcircle", [0x24E5]; "/vdotbelow", [0x1E7F]; "/vecyrillic", [0x0432];
"/veharabic", [0x06A4]; "/vehfinalarabic", [0xFB6B]; "/vehinitialarabic",
[0xFB6C]; "/vehmedialarabic", [0xFB6D]; "/vekatakana", [0x30F9]; "/venus",
[0x2640]; "/verticalbar", [0x007C]; "/verticallineabovecmb", [0x030D];
"/verticallinebelowcmb", [0x0329]; "/verticallinelowmod", [0x02CC];
"/verticallinemod", [0x02C8]; "/vewarmenian", [0x057E]; "/vhook", [0x028B];
"/vikatakana", [0x30F8]; "/viramabengali", [0x09CD]; "/viramadeva", [0x094D];
"/viramagujarati", [0x0ACD]; "/visargabengali", [0x0983]; "/visargadeva",
[0x0903]; "/visargagujarati", [0x0A83]; "/vmonospace", [0xFF56]; "/voarmenian",
[0x0578]; "/voicediterationhiragana", [0x309E]; "/voicediterationkatakana",
[0x30FE]; "/voicedmarkkana", [0x309B]; "/voicedmarkkanahalfwidth", [0xFF9E];
"/vokatakana", [0x30FA]; "/vparen", [0x24B1]; "/vtilde", [0x1E7D]; "/vturned",
[0x028C]; "/vuhiragana", [0x3094]; "/vukatakana", [0x30F4]; "/w", [0x0077];
"/wacute", [0x1E83]; "/waekorean", [0x3159]; "/wahiragana", [0x308F];
"/wakatakana", [0x30EF]; "/wakatakanahalfwidth", [0xFF9C]; "/wakorean",
[0x3158]; "/wasmallhiragana", [0x308E]; "/wasmallkatakana", [0x30EE];
"/wattosquare", [0x3357]; "/wavedash", [0x301C]; "/wavyunderscorevertical",
[0xFE34]; "/wawarabic", [0x0648]; "/wawfinalarabic", [0xFEEE];
"/wawhamzaabovearabic", [0x0624]; "/wawhamzaabovefinalarabic", [0xFE86];
"/wbsquare", [0x33DD]; "/wcircle", [0x24E6]; "/wcircumflex", [0x0175];
"/wdieresis", [0x1E85]; "/wdotaccent", [0x1E87]; "/wdotbelow", [0x1E89];
"/wehiragana", [0x3091]; "/weierstrass", [0x2118]; "/wekatakana", [0x30F1];
"/wekorean", [0x315E]; "/weokorean", [0x315D]; "/wgrave", [0x1E81];
"/whitebullet", [0x25E6]; "/whitecircle", [0x25CB]; "/whitecircleinverse",
[0x25D9]; "/whitecornerbracketleft", [0x300E];
"/whitecornerbracketleftvertical", [0xFE43]; "/whitecornerbracketright",
[0x300F]; "/whitecornerbracketrightvertical", [0xFE44]; "/whitediamond",
[0x25C7]; "/whitediamondcontainingblacksmalldiamond", [0x25C8];
"/whitedownpointingsmalltriangle", [0x25BF]; "/whitedownpointingtriangle",
[0x25BD]; "/whiteleftpointingsmalltriangle", [0x25C3];
"/whiteleftpointingtriangle", [0x25C1]; "/whitelenticularbracketleft", [0x3016];
"/whitelenticularbracketright", [0x3017]; "/whiterightpointingsmalltriangle",
[0x25B9]; "/whiterightpointingtriangle", [0x25B7]; "/whitesmallsquare",
[0x25AB]; "/whitesmilingface", [0x263A]; "/whitesquare", [0x25A1]; "/whitestar",
[0x2606]; "/whitetelephone", [0x260F]; "/whitetortoiseshellbracketleft",
[0x3018]; "/whitetortoiseshellbracketright", [0x3019];
"/whiteuppointingsmalltriangle", [0x25B5]; "/whiteuppointingtriangle", [0x25B3];
"/wihiragana", [0x3090]; "/wikatakana", [0x30F0]; "/wikorean", [0x315F];
"/wmonospace", [0xFF57]; "/wohiragana", [0x3092]; "/wokatakana", [0x30F2];
"/wokatakanahalfwidth", [0xFF66]; "/won", [0x20A9]; "/wonmonospace", [0xFFE6];
"/wowaenthai", [0x0E27]; "/wparen", [0x24B2]; "/wring", [0x1E98]; "/wsuperior",
[0x02B7]; "/wturned", [0x028D]; "/wynn", [0x01BF]; "/x", [0x0078]; "/xabovecmb",
[0x033D]; "/xbopomofo", [0x3112]; "/xcircle", [0x24E7]; "/xdieresis", [0x1E8D];
"/xdotaccent", [0x1E8B]; "/xeharmenian", [0x056D]; "/xi", [0x03BE];
"/xmonospace", [0xFF58]; "/xparen", [0x24B3]; "/xsuperior", [0x02E3]; "/y",
[0x0079]; "/yaadosquare", [0x334E]; "/yabengali", [0x09AF]; "/yacute", [0x00FD];
"/yadeva", [0x092F]; "/yaekorean", [0x3152]; "/yagujarati", [0x0AAF];
"/yagurmukhi", [0x0A2F]; "/yahiragana", [0x3084]; "/yakatakana", [0x30E4];
"/yakatakanahalfwidth", [0xFF94]; "/yakorean", [0x3151]; "/yamakkanthai",
[0x0E4E]; "/yasmallhiragana", [0x3083]; "/yasmallkatakana", [0x30E3];
"/yasmallkatakanahalfwidth", [0xFF6C]; "/yatcyrillic", [0x0463]; "/ycircle",
[0x24E8]; "/ycircumflex", [0x0177]; "/ydieresis", [0x00FF]; "/ydotaccent",
[0x1E8F]; "/ydotbelow", [0x1EF5]; "/yeharabic", [0x064A]; "/yehbarreearabic",
[0x06D2]; "/yehbarreefinalarabic", [0xFBAF]; "/yehfinalarabic", [0xFEF2];
"/yehhamzaabovearabic", [0x0626]; "/yehhamzaabovefinalarabic", [0xFE8A];
"/yehhamzaaboveinitialarabic", [0xFE8B]; "/yehhamzaabovemedialarabic", [0xFE8C];
"/yehinitialarabic", [0xFEF3]; "/yehmedialarabic", [0xFEF4];
"/yehmeeminitialarabic", [0xFCDD]; "/yehmeemisolatedarabic", [0xFC58];
"/yehnoonfinalarabic", [0xFC94]; "/yehthreedotsbelowarabic", [0x06D1];
"/yekorean", [0x3156]; "/yen", [0x00A5]; "/yenmonospace", [0xFFE5];
"/yeokorean", [0x3155]; "/yeorinhieuhkorean", [0x3186]; "/yerahbenyomohebrew",
[0x05AA]; "/yerahbenyomolefthebrew", [0x05AA]; "/yericyrillic", [0x044B];
"/yerudieresiscyrillic", [0x04F9]; "/yesieungkorean", [0x3181];
"/yesieungpansioskorean", [0x3183]; "/yesieungsioskorean", [0x3182];
"/yetivhebrew", [0x059A]; "/ygrave", [0x1EF3]; "/yhook", [0x01B4];
"/yhookabove", [0x1EF7]; "/yiarmenian", [0x0575]; "/yicyrillic", [0x0457];
"/yikorean", [0x3162]; "/yinyang", [0x262F]; "/yiwnarmenian", [0x0582];
"/ymonospace", [0xFF59]; "/yod", [0x05D9]; "/yoddagesh", [0xFB39];
"/yoddageshhebrew", [0xFB39]; "/yodhebrew", [0x05D9]; "/yodyodhebrew", [0x05F2];
"/yodyodpatahhebrew", [0xFB1F]; "/yohiragana", [0x3088]; "/yoikorean", [0x3189];
"/yokatakana", [0x30E8]; "/yokatakanahalfwidth", [0xFF96]; "/yokorean",
[0x315B]; "/yosmallhiragana", [0x3087]; "/yosmallkatakana", [0x30E7];
"/yosmallkatakanahalfwidth", [0xFF6E]; "/yotgreek", [0x03F3]; "/yoyaekorean",
[0x3188]; "/yoyakorean", [0x3187]; "/yoyakthai", [0x0E22]; "/yoyingthai",
[0x0E0D]; "/yparen", [0x24B4]; "/ypogegrammeni", [0x037A];
"/ypogegrammenigreekcmb", [0x0345]; "/yr", [0x01A6]; "/yring", [0x1E99];
"/ysuperior", [0x02B8]; "/ytilde", [0x1EF9]; "/yturned", [0x028E];
"/yuhiragana", [0x3086]; "/yuikorean", [0x318C]; "/yukatakana", [0x30E6];
"/yukatakanahalfwidth", [0xFF95]; "/yukorean", [0x3160]; "/yusbigcyrillic",
[0x046B]; "/yusbigiotifiedcyrillic", [0x046D]; "/yuslittlecyrillic", [0x0467];
"/yuslittleiotifiedcyrillic", [0x0469]; "/yusmallhiragana", [0x3085];
"/yusmallkatakana", [0x30E5]; "/yusmallkatakanahalfwidth", [0xFF6D];
"/yuyekorean", [0x318B]; "/yuyeokorean", [0x318A]; "/yyabengali", [0x09DF];
"/yyadeva", [0x095F]; "/z", [0x007A]; "/zaarmenian", [0x0566]; "/zacute",
[0x017A]; "/zadeva", [0x095B]; "/zagurmukhi", [0x0A5B]; "/zaharabic", [0x0638];
"/zahfinalarabic", [0xFEC6]; "/zahinitialarabic", [0xFEC7]; "/zahiragana",
[0x3056]; "/zahmedialarabic", [0xFEC8]; "/zainarabic", [0x0632];
"/zainfinalarabic", [0xFEB0]; "/zakatakana", [0x30B6]; "/zaqefgadolhebrew",
[0x0595]; "/zaqefqatanhebrew", [0x0594]; "/zarqahebrew", [0x0598]; "/zayin",
[0x05D6]; "/zayindagesh", [0xFB36]; "/zayindageshhebrew", [0xFB36];
"/zayinhebrew", [0x05D6]; "/zbopomofo", [0x3117]; "/zcaron", [0x017E];
"/zcircle", [0x24E9]; "/zcircumflex", [0x1E91]; "/zcurl", [0x0291]; "/zdot",
[0x017C]; "/zdotaccent", [0x017C]; "/zdotbelow", [0x1E93]; "/zecyrillic",
[0x0437]; "/zedescendercyrillic", [0x0499]; "/zedieresiscyrillic", [0x04DF];
"/zehiragana", [0x305C]; "/zekatakana", [0x30BC]; "/zero", [0x0030];
"/zeroarabic", [0x0660]; "/zerobengali", [0x09E6]; "/zerodeva", [0x0966];
"/zerogujarati", [0x0AE6]; "/zerogurmukhi", [0x0A66]; "/zerohackarabic",
[0x0660]; "/zeroinferior", [0x2080]; "/zeromonospace", [0xFF10];
"/zerooldstyle", [0xF730]; "/zeropersian", [0x06F0]; "/zerosuperior", [0x2070];
"/zerothai", [0x0E50]; "/zerowidthjoiner", [0xFEFF]; "/zerowidthnonjoiner",
[0x200C]; "/zerowidthspace", [0x200B]; "/zeta", [0x03B6]; "/zhbopomofo",
[0x3113]; "/zhearmenian", [0x056A]; "/zhebrevecyrillic", [0x04C2];
"/zhecyrillic", [0x0436]; "/zhedescendercyrillic", [0x0497];
"/zhedieresiscyrillic", [0x04DD]; "/zihiragana", [0x3058]; "/zikatakana",
[0x30B8]; "/zinorhebrew", [0x05AE]; "/zlinebelow", [0x1E95]; "/zmonospace",
[0xFF5A]; "/zohiragana", [0x305E]; "/zokatakana", [0x30BE]; "/zparen", [0x24B5];
"/zretroflexhook", [0x0290]; "/zstroke", [0x01B6]; "/zuhiragana", [0x305A];
"/zukatakana", [0x30BA]]


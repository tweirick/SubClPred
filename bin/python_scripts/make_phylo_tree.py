from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle

def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=30)
        faces.add_face_to_node(N, node, 0, position="aligned")

def get_example_tree():

    # Set dashed blue lines in all leaves
    #nst1 = NodeStyle()
    #nst1["bgcolor"] = "LightSteelBlue"
    #nst2 = NodeStyle()
    #nst2["bgcolor"] = "Moccasin"
    #nst3 = NodeStyle()
    #nst3["bgcolor"] = "DarkSeaGreen"
    #nst4 = NodeStyle()
    #nst4["bgcolor"] = "Khaki"





    t = Tree("""( ( ( ( sp|Q02075|LAC2_THACU:0.276276 , ( sp|Q02081|LAC4_THACU:0.247646 , tr|I3PL63|I3PL63_FLAVE:0.247646 ):0.0286296 ):0.0355952 , ( ( ( tr|Q4VY49|Q4VY49_9AGAR:0.0105566 , ( tr|Q7Z8S4|Q7Z8S4_PLESA:0.000959693 , tr|Q2VT19|Q2VT19_PLEPU:0.000959693 ):0.00959693 ):0.155549 , ( sp|Q12541|LAC1_AGABI:0.0346154 , sp|Q12542|LAC2_AGABI:0.0346154 ):0.13149 ):0.136677 , ( ( ( tr|Q2AAD1|Q2AAD1_FLAVE:0.179264 , tr|G8A520|G8A520_FLAVE:0.179264 ):0.0840652 , ( ( tr|Q69FW9|Q69FW9_9AGAR:0.162295 , ( tr|Q69FW8|Q69FW8_9AGAR:0.155018 , ( tr|Q69FX1|Q69FX1_9AGAR:0.115494 , ( tr|Q6X938|Q6X938_9AGAR:0.0935039 , tr|Q69FW7|Q69FW7_9AGAR:0.0935039 ):0.0219899 ):0.0395246 ):0.00727676 ):0.0985383 , ( tr|Q9P8B9|Q9P8B9_COPCO:0.259043 , ( tr|Q9Y780|Q9Y780_COPCI:0.24834 , ( tr|O74171|O74171_SCHCO:0.212551 , ( tr|Q69FX0|Q69FX0_9AGAR:0.206553 , ( ( ( ( tr|B5MAF5|B5MAF5_PLEOS:0.0422932 , tr|Q7Z8S5|Q7Z8S5_PLESA:0.0422932 ):0.0906258 , ( tr|D2KZ02|D2KZ02_9AGAR:0.135248 , ( ( tr|Q6A1A1|Q6A1A1_9AGAR:0.0338983 , ( ( sp|Q12739|LAC2_PLEOS:0.0145403 , ( tr|Q9UVY4|Q9UVY4_PLEOS:0.00375235 , tr|Q6STF0|Q6STF0_PLEOS:0.00375235 ):0.010788 ):0.0122454 , ( tr|Q7Z8S3|Q7Z8S3_PLESA:0.0056391 , tr|Q2VT18|Q2VT18_PLEPU:0.0056391 ):0.0211466 ):0.00711259 ):0.0167037 , ( ( tr|D3YJ58|D3YJ58_PLEER:0.00376648 , tr|B0JDP8|B0JDP8_PLEER:0.00376648 ):0.0398573 , ( tr|Q7Z8S6|Q7Z8S6_PLESA:0.0344991 , ( sp|Q12729|LAC1_PLEOS:0.00283554 , tr|Q6RYA4|Q6RYA4_PLEOS:0.00283554 ):0.0316635 ):0.00912472 ):0.00697827 ):0.0846458 ):-0.0023288 ):0.0408624 , ( tr|Q7Z8S2|Q7Z8S2_PLESA:0.141429 , ( tr|D2KZ06|D2KZ06_9AGAR:0.00471698 , tr|D2KZ07|D2KZ07_9AGAR:0.00471698 ):0.136712 ):0.0323528 ):0.038193 , ( tr|C0JRG5|C0JRG5_LACBI:0.211952 , ( tr|G8A545|G8A545_FLAVE:0.204146 , ( tr|H9BT70|H9BT70_9AGAR:0.186331 , ( ( tr|H9BT71|H9BT71_9AGAR:0.125486 , tr|A6N8N5|A6N8N5_PHONA:0.125486 ):0.0417706 , ( tr|A8W7J6|A8W7J6_9AGAR:0.159596 , ( tr|H9C4A3|H9C4A3_COPCM:0.176783 , ( tr|H9C4A2|H9C4A2_COPCM:0.178283 , ( tr|C0JRH0|C0JRH0_LACBI:0.170859 , ( tr|C0JRG6|C0JRG6_LACBI:0.109373 , ( tr|C0JRG4|C0JRG4_LACBI:0.0822594 , ( tr|C0JRG3|C0JRG3_LACBI:0.0735587 , tr|C0JRG2|C0JRG2_LACBI:0.0735587 ):0.00870072 ):0.0271135 ):0.0614862 ):0.00742415 ):-0.00150046 ):-0.017187 ):0.00766109 ):0.0190739 ):0.0178153 ):0.00780588 ):2.22921e-05 ):-0.0054214 ):0.00599766 ):0.0357894 ):0.0107031 ):0.00179023 ):0.00249535 ):0.0335869 , ( tr|G8A555|G8A555_FLAVE:0.284044 , ( ( tr|C5NN27|C5NN27_LENED:0.217025 , ( ( tr|Q8X1W2|Q8X1W2_LENED:0.00638686 , tr|Q68LM0|Q68LM0_LENED:0.00638686 ):0.179258 , ( tr|Q8WZH9|Q8WZH9_LENED:0.0181421 , ( tr|Q68LM3|Q68LM3_LENED:0.0146067 , ( tr|Q8WZI0|Q8WZI0_LENED:0.00891799 , ( tr|Q68LM4|Q68LM4_LENED:0.00686284 , ( tr|Q68LM2|Q68LM2_LENED:0.00506073 , tr|Q68LM5|Q68LM5_LENED:0.00506073 ):0.00180211 ):0.00205515 ):0.00568875 ):0.00353539 ):0.167503 ):0.0313805 ):0.0430473 , ( ( ( tr|B5MAF4|B5MAF4_PLEOS:0.113985 , tr|D2KZ05|D2KZ05_9AGAR:0.113985 ):0.115398 , ( tr|D2KZ04|D2KZ04_9AGAR:0.217866 , ( ( tr|D2KZ01|D2KZ01_9AGAR:0.101799 , ( tr|D2KZ00|D2KZ00_9AGAR:0.0047259 , tr|D2KZ03|D2KZ03_9AGAR:0.0047259 ):0.0970733 ):0.0964945 , ( tr|D2KYZ9|D2KYZ9_9AGAR:0.154986 , ( ( tr|Q308C0|Q308C0_PLEER:0.00140713 , ( tr|B6V331|B6V331_PLEER:0.000938086 , tr|B0JDP9|B0JDP9_PLEER:0.000938086 ):0.000469043 ):0.0214927 , ( tr|O60199|O60199_PLEOS:0.00972763 , tr|D4AIA5|D4AIA5_PLEOS:0.00972763 ):0.0131722 ):0.132086 ):0.043308 ):0.0195725 ):0.0115165 ):0.0137794 , ( tr|G8A529|G8A529_FLAVE:0.247393 , ( ( tr|C0JRG7|C0JRG7_LACBI:0.190554 , ( tr|C0JRG9|C0JRG9_LACBI:0.112186 , tr|C0JRG8|C0JRG8_LACBI:0.112186 ):0.0783681 ):0.0376743 , ( tr|Q6H9H7|Q6H9H7_9APHY:0.223384 , ( tr|Q9P8G4|Q9P8G4_9APHY:0.204316 , ( ( tr|I1SB14|I1SB14_9APHY:0.10303 , tr|I1VE66|I1VE66_9APHY:0.10303 ):0.0900954 , ( tr|I1VE67|I1VE67_9APHY:0.184866 , ( tr|B8YQ97|B8YQ97_CERUI:0.129412 , tr|Q1EPM3|Q1EPM3_9APHY:0.129412 ):0.0554546 ):0.00825934 ):0.0111907 ):0.0190677 ):0.00484392 ):0.0191652 ):-0.00423114 ):0.0169105 ):0.0239713 ):0.0128717 ):0.0058676 ):0.00908756 ):0.353624 , ( ( tr|I1W1V8|I1W1V8_9APHY:0.15251 , ( tr|O59944|O59944_CERSU:0.0693642 , tr|I1W1V7|I1W1V7_9APHY:0.0693642 ):0.0831455 ):0.0230437 , ( ( sp|Q01679|LAC1_PHLRA:0.150968 , ( ( tr|E7BLQ8|E7BLQ8_9APHY:0.0705996 , tr|E7BLQ9|E7BLQ9_9APHY:0.0705996 ):0.0736765 , ( tr|Q6UNT7|Q6UNT7_9APHY:0.128627 , tr|I1VE65|I1VE65_9APHY:0.128627 ):0.0156494 ):0.00669187 ):0.023287 , ( ( tr|G9I8W7|G9I8W7_LENED:0.00361969 , ( tr|Q68LM1|Q68LM1_LENED:0.00337838 , ( tr|Q8WZG3|Q8WZG3_LENED:0.000965251 , tr|G9I8W6|G9I8W6_LENED:0.000965251 ):0.00241313 ):0.000241313 ):0.170713 , ( ( tr|E7BLR0|E7BLR0_9APHY:0.0910853 , tr|Q2HWK1|Q2HWK1_9APHY:0.0910853 ):0.0957099 , ( ( tr|Q50JG4|Q50JG4_TRAVE:0.13125 , ( tr|Q5I7J0|Q5I7J0_9APHY:0.0780952 , ( sp|Q12717|LAC5_TRAVE:0.0104364 , tr|Q50JG3|Q50JG3_TRAVE:0.0104364 ):0.0676588 ):0.0531548 ):0.0354404 , ( ( ( tr|Q8NID5|Q8NID5_9APHY:0.0284091 , tr|C5IXN8|C5IXN8_9APHY:0.0284091 ):0.129168 , ( tr|C1KDZ5|C1KDZ5_9APHY:0.122104 , ( tr|C1KDZ8|C1KDZ8_9APHY:0.0780347 , ( tr|C1KDZ6|C1KDZ6_9APHY:0.0634615 , ( tr|Q6TH77|Q6TH77_9APHY:0.0400763 , tr|C1KDZ7|C1KDZ7_9APHY:0.0400763 ):0.0233852 ):0.0145731 ):0.0440696 ):0.0354729 ):0.0146112 , ( ( ( tr|Q6RYA5|Q6RYA5_FLAVE:0.00191939 , tr|G0WM60|G0WM60_9APHY:0.00191939 ):0.109165 , ( tr|B5G556|B5G556_9APHY:0.0494242 , ( tr|B5G554|B5G554_GANTS:0.0143954 , tr|K9R5B2|K9R5B2_GANLU:0.0143954 ):0.0350288 ):0.0616603 ):0.042677 , ( ( ( tr|Q6R5P8|Q6R5P8_9APHY:0.0212355 , tr|C1JCL7|C1JCL7_9APHY:0.0212355 ):0.135561 , ( ( tr|G9M4T7|G9M4T7_TRAVE:0.0639313 , tr|Q716A0|Q716A0_9APHY:0.0639313 ):0.067028 , ( ( tr|M1GME7|M1GME7_9APHY:0.0408607 , ( ( tr|Q1W6B1|Q1W6B1_9APHY:0.00580271 , tr|C6G7V1|C6G7V1_9APHY:0.00580271 ):0.0140232 , ( tr|Q9HDQ0|Q9HDQ0_9APHY:0.0145068 , tr|Q9UVQ5|Q9UVQ5_9APHY:0.0145068 ):0.00531915 ):0.0210348 ):0.093711 , ( ( tr|Q9UVQ2|Q9UVQ2_PYCCI:0.0458494 , ( sp|O59896|LAC1_PYCCI:0.0236486 , ( tr|C9WKP8|C9WKP8_9APHY:0.011583 , tr|Q96TR6|Q96TR6_PYCCO:0.011583 ):0.0120656 ):0.0222008 ):0.072212 , ( tr|Q50JG6|Q50JG6_TRAVE:0.123014 , ( ( sp|D0VWU3|LAC1_TRAMX:0.0360721 , ( tr|B2L9C1|B2L9C1_TRAHI:0.0184466 , tr|Q8TFL8|Q8TFL8_TRAHI:0.0184466 ):0.0176255 ):0.0409032 , ( tr|Q716A1|Q716A1_9APHY:0.0832031 , ( ( sp|Q02497|LAC1_TRAHI:0.0448317 , ( tr|Q8TG94|Q8TG94_TRAPU:0.0175481 , ( tr|O94222|O94222_TRAVE:0.00865385 , ( tr|Q96UT7|Q96UT7_TRAVE:0.00576923 , tr|F4ZCH1|F4ZCH1_9APHY:0.00576923 ):0.00288462 ):0.00889423 ):0.0272837 ):0.0409405 , ( tr|G4XIH4|G4XIH4_9APHY:0.0754808 , ( tr|G4XU43|G4XU43_9APHY:0.0471154 , tr|D9J137|D9J137_9APHY:0.0471154 ):0.0283654 ):0.0102915 ):-0.00256911 ):-0.00622776 ):0.0460388 ):-0.00495272 ):0.0165103 ):-0.0036124 ):0.0258374 ):-0.0237401 , ( ( ( tr|Q716A3|Q716A3_9APHY:0.00192308 , tr|Q716A2|Q716A2_9APHY:0.00192308 ):0.0716161 , ( tr|F2VPT7|F2VPT7_9APHY:0.0276975 , ( tr|Q96UK8|Q96UK8_TRAVE:0.0134875 , ( tr|B8Y3J5|B8Y3J5_TRAVE:0.0115607 , ( tr|Q50JG5|Q50JG5_TRAVE:0.00963391 , ( sp|Q12718|LAC2_TRAVE:0.00578035 , tr|Q5IR80|Q5IR80_TRAVE:0.00578035 ):0.00385357 ):0.00192678 ):0.00192678 ):0.01421 ):0.0458417 ):0.0474636 , ( ( tr|I3NL60|I3NL60_9APHY:0.0703123 , ( tr|E1U755|E1U755_9APHY:0 , tr|E1U754|E1U754_9APHY:0 ):0.0703123 ):0.0413365 , ( tr|A3F8Z8|A3F8Z8_9APHY:0.100962 , ( tr|Q9HG17|Q9HG17_GANLU:0.0711538 , tr|B5G552|B5G552_GANLU:0.0711538 ):0.0298077 ):0.0106873 ):0.009354 ):0.0120539 ):0.0207048 ):0.018427 ):-0.00549796 ):0.0201047 ):-0.0124622 ):-7.80523e-05 ):0.00129837 ):0.353624 ):0.82517 , ( ( ( tr|Q8W0V5|Q8W0V5_LOLPR:0.260525 , ( sp|Q339K6|LAC15_ORYSJ:0.242491 , ( sp|Q53LU4|LAC18_ORYSJ:0.0930413 , ( ( sp|A2Y9C5|LAC19_ORYSI:0.00169492 , sp|Q2R0L2|LAC19_ORYSJ:0.00169492 ):0.0610206 , ( sp|A2Y9C2|LAC20_ORYSI:0.00862069 , sp|Q2R0L0|LAC20_ORYSJ:0.00862069 ):0.0540948 ):0.0303258 ):0.149449 ):0.0180341 ):0.067436 , ( ( ( ( sp|Q9LFD2|LAC8_ARATH:0.0590753 , sp|Q9LFD1|LAC9_ARATH:0.0590753 ):0.174155 , ( tr|Q2LD62|Q2LD62_PEA:0.224574 , ( sp|Q9SR40|LAC7_ARATH:0.187831 , tr|B1PXG7|B1PXG7_SOLLC:0.187831 ):0.0367428 ):0.00865725 ):0.0485124 , ( ( sp|Q69L99|LAC14_ORYSJ:0.249141 , ( sp|Q5N7A3|LAC6_ORYSJ:0.211063 , ( sp|Q2QUN2|LAC24_ORYSJ:0.159445 , sp|Q0IP28|LAC25_ORYSJ:0.159445 ):0.0516174 ):0.0380785 ):0.0228225 , ( sp|Q0JHP8|LAC8_ORYSJ:0.251579 , ( tr|Q4VJ26|Q4VJ26_MAIZE:0.165633 , ( sp|Q5N7B4|LAC7_ORYSJ:0.132379 , tr|Q8W0V4|Q8W0V4_LOLPR:0.132379 ):0.0332539 ):0.0859463 ):0.0203843 ):0.00977939 ):0.026731 , ( sp|Q9ZPY2|LAC6_ARATH:0.303579 , ( ( sp|Q9LMS3|LAC1_ARATH:0.269856 , ( ( ( tr|Q2PAJ2|Q2PAJ2_MAIZE:0.159796 , ( sp|Q5N9W4|LAC5_ORYSJ:0.0941499 , tr|Q2PAI9|Q2PAI9_MAIZE:0.0941499 ):0.0656463 ):0.0626486 , ( ( tr|Q9AUI6|Q9AUI6_PINTA:0.153325 , ( tr|Q9AUI1|Q9AUI1_PINTA:0.0865052 , tr|Q9AUI2|Q9AUI2_PINTA:0.0865052 ):0.0668201 ):0.037825 , ( tr|M5AP95|M5AP95_CHAOB:0.194695 , ( ( ( tr|Q8W0V6|Q8W0V6_LOLPR:0.120675 , ( sp|Q5N9X2|LAC4_ORYSJ:0.0906736 , tr|Q2PAJ0|Q2PAJ0_MAIZE:0.0906736 ):0.0300012 ):0.0158563 , ( sp|Q0DHL2|LAC12_ORYSJ:0.108885 , tr|Q8S2A8|Q8S2A8_ORYSJ:0.108885 ):0.027646 ):0.0473467 , ( tr|Q9AUI4|Q9AUI4_PINTA:0.185505 , ( sp|O81081|LAC2_ARATH:0.180085 , ( sp|Q10ND7|LAC10_ORYSJ:0.163175 , ( sp|Q9FJD5|LAC17_ARATH:0.151703 , ( ( tr|Q9ZQW2|Q9ZQW2_POPTR:0.00431034 , tr|B9HBS9|B9HBS9_POPTR:0.00431034 ):0.145261 , ( tr|O24041|O24041_LIRTU:0.118421 , ( tr|O24043|O24043_LIRTU:0.0653846 , ( tr|O24044|O24044_LIRTU:0.0444444 , tr|O24042|O24042_LIRTU:0.0444444 ):0.0209402 ):0.0530364 ):0.0311498 ):0.00213183 ):0.0114721 ):0.0169104 ):0.00541964 ):-0.0016271 ):0.0108176 ):-0.00354503 ):0.0312944 ):0.0224021 , ( ( ( tr|B9HHK7|B9HHK7_POPTR:0.0551802 , ( tr|B9HT10|B9HT10_POPTR:0.0130631 , ( tr|Q9FSC9|Q9FSC9_POPTR:0.000900901 , tr|Q9ZQW3|Q9ZQW3_POPTR:0.000900901 ):0.0121622 ):0.0421171 ):0.13467 , ( sp|Q1PDH6|LAC16_ARATH:0.178506 , ( sp|Q0IQU1|LAC22_ORYSJ:0.175137 , ( sp|Q6ID18|LAC10_ARATH:0.13475 , ( sp|O80434|LAC4_ARATH:0.117594 , ( tr|B9IG56|B9IG56_POPTR:0.0816876 , tr|P93366|P93366_TOBAC:0.0816876 ):0.0359066 ):0.0171562 ):0.0403865 ):0.003369 ):0.0113439 ):0.0178583 , ( ( tr|Q9AUH9|Q9AUH9_PINTA:0.142793 , ( tr|Q9AUI0|Q9AUI0_PINTA:0.132432 , tr|M5AN30|M5AN30_CHAOB:0.132432 ):0.0103604 ):0.0544501 , ( sp|Q8VZA1|LAC11_ARATH:0.155296 , sp|Q8RYM9|LAC2_ORYSJ:0.155296 ):0.0419467 ):0.0104653 ):0.0371387 ):0.025009 ):0.0151808 , ( ( ( tr|K4P1L9|K4P1L9_PICAB:0.00395431 , ( tr|K4P1P7|K4P1P7_PICAB:0.00263158 , tr|K4NZ22|K4NZ22_PICAB:0.00263158 ):0.00132273 ):0.229173 , ( ( sp|Q9LYQ2|LAC13_ARATH:0.139719 , sp|Q56YT0|LAC3_ARATH:0.139719 ):0.0684609 , ( ( ( tr|K4P1M3|K4P1M3_PICAB:0.0008726 , tr|K4PCR3|K4PCR3_PICAB:0.0008726 ):0.146564 , ( ( tr|K4NZE7|K4NZE7_PICAB:0.0701754 , ( tr|Q9AUI3|Q9AUI3_PINTA:0.0649123 , tr|K4PCQ7|K4PCQ7_PICAB:0.0649123 ):0.00526316 ):0.0224918 , ( tr|K4NZF4|K4NZF4_PICAB:0.09277 , ( tr|Q9AUI5|Q9AUI5_PINTA:0.0199653 , tr|F4MKL7|F4MKL7_PINPS:0.0199653 ):0.0728048 ):-0.000102796 ):0.0547689 ):0.0489896 , ( sp|Q2QYS3|LAC23_ORYSJ:0.180526 , ( ( sp|Q941X2|LAC3_ORYSJ:0.0767196 , tr|C0P5Q0|C0P5Q0_MAIZE:0.0767196 ):0.0907467 , ( ( tr|Q8L4Y3|Q8L4Y3_SOYBN:0.113676 , ( tr|Q9ZP47|Q9ZP47_POPTR:0.0121951 , tr|B9HHV7|B9HHV7_POPTR:0.0121951 ):0.101481 ):0.0117478 , ( sp|Q9FLB5|LAC12_ARATH:0.111504 , sp|Q9SIY8|LAC5_ARATH:0.111504 ):0.0139193 ):0.0420425 ):0.0130598 ):0.0158997 ):0.0117539 ):0.0249472 ):0.0541989 , ( ( tr|Q6TDS6|Q6TDS6_GOSAR:0.0344523 , tr|G8XQW0|G8XQW0_GOSHI:0.0344523 ):0.247136 , ( sp|Q9FY79|LAC14_ARATH:0.268169 , ( ( tr|C9E6Q3|C9E6Q3_9ROSI:0.201582 , ( tr|Q38757|Q38757_ACEPS:0.113274 , tr|B2CMA7|B2CMA7_LITCN:0.113274 ):0.0883072 ):0.0578095 , ( ( sp|Q84J37|LAC15_ARATH:0.0970228 , ( tr|G0WXI9|G0WXI9_BRANA:0.0786449 , ( tr|G0WXI5|G0WXI5_BRANA:0.0268336 , tr|G8Z904|G8Z904_BRANA:0.0268336 ):0.0518113 ):0.0183779 ):0.157435 , ( tr|B2M153|B2M153_ROSHC:0.238438 , ( tr|Q2PAJ1|Q2PAJ1_MAIZE:0.191451 , ( sp|Q2QZ80|LAC21_ORYSJ:0.125216 , sp|Q6Z8L2|LAC9_ORYSJ:0.125216 ):0.0662355 ):0.0469867 ):0.0160202 ):0.00493276 ):0.00877827 ):0.0134194 ):0.00573713 ):-0.00228906 ):0.018542 ):0.00489542 ):0.0194865 ):0.776588 , ( ( ( tr|J9PBQ8|J9PBQ8_9PSEU:0.200717 , tr|J9PBR2|J9PBR2_STRVD:0.200717 ):0.233755 , ( ( sp|P17489|LAC1_EMENI:0.392387 , ( tr|Q6STE9|Q6STE9_9HOMO:0.0634648 , tr|D3K4I1|D3K4I1_9HOMO:0.0634648 ):0.328922 ):0.0346513 , ( tr|Q72HW2|Q72HW2_THET2:0.435706 , ( tr|Q5N7B3|Q5N7B3_ORYSJ:0.428232 , ( ( tr|Q9VX11|Q9VX11_DROME:0.360428 , ( ( tr|Q8I8Y2|Q8I8Y2_ANOGA:0.29623 , ( ( tr|D5MRE2|D5MRE2_NEPCI:0.241084 , tr|D5MRE1|D5MRE1_NEPCI:0.241084 ):0.0381306 , ( tr|Q8I8Y1|Q8I8Y1_MANSE:0.266453 , ( tr|Q49I37|Q49I37_TRICA:0.239598 , tr|M4GPQ6|M4GPQ6_9NEOP:0.239598 ):0.0268549 ):0.0127617 ):0.0170157 ):0.0490155 , ( ( tr|Q4U3X4|Q4U3X4_AEDAE:0.207262 , tr|A5YVV0|A5YVV0_ANOGA:0.207262 ):0.129277 , ( ( ( tr|D0E8H3|D0E8H3_RETFL:0.000959693 , ( tr|D0E8H1|D0E8H1_RETFL:0.0011592 , ( tr|D0E8H7|D0E8H7_RETFL:0.000772798 , tr|D0E8H0|D0E8H0_RETFL:0.000772798 ):0.000386399 ):-0.000199503 ):0.013182 , ( tr|D0E8H2|D0E8H2_RETFL:0.00237083 , ( ( tr|D0E8H4|D0E8H4_RETFL:0.0015456 , tr|D0E8H5|D0E8H5_RETFL:0.0015456 ):6.99139e-05 , tr|D0E8H6|D0E8H6_RETFL:0.00161551 ):0.000755319 ):0.0117709 ):0.300337 , ( tr|Q8WPD1|Q8WPD1_PIMHY:0.289625 , ( ( tr|D5MRE3|D5MRE3_NEPCI:0.0969784 , ( tr|E9RH10|E9RH10_9HEMI:0.057931 , tr|E9RH11|E9RH11_9HEMI:0.057931 ):0.0390473 ):0.0513012 , ( ( ( tr|D5MRS1|D5MRS1_9NEOP:0.0129482 , ( tr|D4AH59|D4AH59_9NEOP:0.00464807 , tr|D5MRR9|D5MRR9_9NEOP:0.00464807 ):0.00830013 ):0.0639318 , ( tr|B4F7L6|B4F7L6_BOMMO:0.0668143 , ( tr|Q8I8Y0|Q8I8Y0_MANSE:0.0319079 , ( tr|A7XQR9|A7XQR9_BOMMO:0.00262123 , tr|B5BR55|B5BR55_BOMMO:0.00262123 ):0.0292867 ):0.0349065 ):0.0100656 ):0.0554303 , ( ( ( tr|A1Z6F6|A1Z6F6_DROME:0.00600801 , tr|Q7JQF6|Q7JQF6_DROME:0.00600801 ):0.0951285 , ( ( tr|Q58IU3|Q58IU3_ANOGA:0.0332005 , tr|Q58IU2|Q58IU2_ANOGA:0.0332005 ):0.0437187 , ( tr|B5B2D0|B5B2D0_CULPI:0.0475234 , tr|Q4ZGM4|Q4ZGM4_AEDAE:0.0475234 ):0.0293958 ):0.0242173 ):0.0298379 , ( tr|I0IV86|I0IV86_GRYBI:0.11955 , ( tr|A7XQS2|A7XQS2_MONAT:0.107354 , ( tr|Q49I41|Q49I41_TRICA:0.0420757 , tr|Q49I40|Q49I40_TRICA:0.0420757 ):0.0652785 ):0.0121957 ):0.0114244 ):0.00133589 ):0.0159694 ):0.141346 ):0.0248528 ):0.0220606 ):0.00870708 ):0.0151823 ):0.0673417 , ( tr|I3W7E6|I3W7E6_BETVU:0.418472 , ( ( tr|Q6E124|Q6E124_9HOMO:0.401472 , ( ( sp|Q12570|LAC1_BOTFU:0.238414 , ( tr|H8ZRU2|H8ZRU2_9HELO:0.0362069 , sp|Q96WM9|LAC2_BOTFU:0.0362069 ):0.202207 ):0.0919844 , ( ( sp|P06811|LAC1_NEUCR:0.198725 , ( sp|P78722|LAC2_PODAS:0.167816 , ( sp|Q70KY3|LAC1_MELAO:0.140728 , tr|F6N9E7|F6N9E7_9PEZI:0.140728 ):0.0270873 ):0.0309089 ):0.121478 , ( tr|B3TLA6|B3TLA6_9PLEO:0.303333 , tr|I0AVQ6|I0AVQ6_COCLU:0.303333 ):0.0168693 ):0.0101953 ):0.0710741 ):-0.0014559 , ( tr|Q09HV5|Q09HV5_CRYNV:0.407914 , ( ( sp|P56193|LAC1_THACU:0.0664336 , sp|Q02079|LAC3_THACU:0.0664336 ):0.281474 , ( tr|G8A560|G8A560_FLAVE:0.165996 , tr|G8A542|G8A542_FLAVE:0.165996 ):0.181911 ):0.0600065 ):-0.00789782 ):0.0184559 ):0.00929803 ):0.000462234 ):0.00747371 ):-0.00866809 ):0.00743428 ):0.0108743 , ( ( tr|F8V189|F8V189_POPTR:0.332461 , tr|F8V190|F8V190_POPTR:0.332461 ):0.10854 , ( sp|Q99US8|Y1187_STAAM:0.410059 , ( sp|P33644|YFIH_ECOLI:0.398148 , ( sp|Q8BZT9|LACC1_MOUSE:0.115988 , ( tr|K9IXR5|K9IXR5_DESRO:0.0784884 , ( sp|Q8IV20|LACC1_HUMAN:0.0197674 , tr|F6UMP1|F6UMP1_MACMU:0.0197674 ):0.0587209 ):0.0375 ):0.28216 ):0.0119105 ):0.0309426 ):0.00434521 ):0.776588 ):0.402206 ) ;""")

    for n in t.traverse():
        n.dist = 0
    
    #n1 = t.get_common_ancestor("a1", "a2", "a3")
    #n1.set_style(nst1)
    #n2 = t.get_common_ancestor("b1", "b2", "b3", "b4")
    #n2.set_style(nst2)
    #n3 = t.get_common_ancestor("c1", "c2", "c3")
    #n3.set_style(nst3)
    #n4 = t.get_common_ancestor("b3", "b4")
    #n4.set_style(nst4)
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False

    ts.mode = "c"
    ts.root_opening_factor = 1
    return t, ts
    
if __name__ == "__main__":
    t, ts = get_example_tree()
    #t.render("node_background.png", w=400, tree_style=ts)
    t.show(tree_style=ts)

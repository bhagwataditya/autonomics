# scalar
x <- 'Reporter intensity 1 WT(126).KD(127).R1'
dequantify(x)
demultiplex(dequantify(x))

# vector
x <- c('Reporter intensity 1 WT(126).KD(127).R1', 
       'Reporter intensity 2 WT(126).KD(127).R1')
dequantify(x)
demultiplex(dequantify(x))

# blank channels
x <- c('Reporter intensity 1 WT(126).KD(127).R1', 
       'Reporter intensity 2 WT(126).KD(127).R1', 
       'Reporter intensity 3 WT(126).KD(127).R1')
dequantify(x)
demultiplex(dequantify(x))

# replicate embedded in sample name
x <- c('Reporter intensity 1 WT.R1(126).KD.R1(127)', 
       'Reporter intensity 2 WT.R1(126).KD.R1(127)', 
       'Reporter intensity 3 WT.R1(126).KD.R1(127)')
dequantify(x)
demultiplex(dequantify(x))


# 
x <- c(
    'Reporter intensity corrected 1 K3941(126).K3942(127).K3943(128).K3944(129).K3945(130).K3946(131).R1', 
    'Reporter intensity corrected 2 K3941(126).K3942(127).K3943(128).K3944(129).K3945(130).K3946(131).R1', 
    'Reporter intensity corrected 3 K3941(126).K3942(127).K3943(128).K3944(129).K3945(130).K3946(131).R1', 
    'Reporter intensity corrected 4 K3941(126).K3942(127).K3943(128).K3944(129).K3945(130).K3946(131).R1', 
    'Reporter intensity corrected 5 K3941(126).K3942(127).K3943(128).K3944(129).K3945(130).K3946(131).R1', 
    'Reporter intensity corrected 6 K3941(126).K3942(127).K3943(128).K3944(129).K3945(130).K3946(131).R1', 
    
    'Reporter intensity corrected 1 K3941(126).K3942(127).K3947(128).K3948(129).K3949(130).R2', 
    'Reporter intensity corrected 2 K3941(126).K3942(127).K3947(128).K3948(129).K3949(130).R2', 
    'Reporter intensity corrected 3 K3941(126).K3942(127).K3947(128).K3948(129).K3949(130).R2', 
    'Reporter intensity corrected 4 K3941(126).K3942(127).K3947(128).K3948(129).K3949(130).R2', 
    'Reporter intensity corrected 5 K3941(126).K3942(127).K3947(128).K3948(129).K3949(130).R2', 
    'Reporter intensity corrected 6 K3941(126).K3942(127).K3947(128).K3948(129).K3949(130).R2', 
    
    'Reporter intensity corrected 1 K3951(126).K3952(127).K3953(128).K3954(129).K3955(130).K3956(131).R3', 
    'Reporter intensity corrected 2 K3951(126).K3952(127).K3953(128).K3954(129).K3955(130).K3956(131).R3', 
    'Reporter intensity corrected 3 K3951(126).K3952(127).K3953(128).K3954(129).K3955(130).K3956(131).R3', 
    'Reporter intensity corrected 4 K3951(126).K3952(127).K3953(128).K3954(129).K3955(130).K3956(131).R3', 
    'Reporter intensity corrected 5 K3951(126).K3952(127).K3953(128).K3954(129).K3955(130).K3956(131).R3', 
    'Reporter intensity corrected 6 K3951(126).K3952(127).K3953(128).K3954(129).K3955(130).K3956(131).R3', 
    
    'Reporter intensity corrected 1 K3951(126).K3952(127).K3957(128).K3958(129).K3959(130).R4', 
    'Reporter intensity corrected 2 K3951(126).K3952(127).K3957(128).K3958(129).K3959(130).R4', 
    'Reporter intensity corrected 3 K3951(126).K3952(127).K3957(128).K3958(129).K3959(130).R4', 
    'Reporter intensity corrected 4 K3951(126).K3952(127).K3957(128).K3958(129).K3959(130).R4', 
    'Reporter intensity corrected 5 K3951(126).K3952(127).K3957(128).K3958(129).K3959(130).R4', 
    'Reporter intensity corrected 6 K3951(126).K3952(127).K3957(128).K3958(129).K3959(130).R4', 
    
    'Reporter intensity corrected 1 K3971(126).K3972(127).K39711(128).BLANK(129).K3975(130).K39710(131).R7', 
    'Reporter intensity corrected 2 K3971(126).K3972(127).K39711(128).BLANK(129).K3975(130).K39710(131).R7', 
    'Reporter intensity corrected 3 K3971(126).K3972(127).K39711(128).K3975(130).K39710(131).R7', 
    'Reporter intensity corrected 4 K3971(126).K3972(127).K39711(128).K3975(130).K39710(131).R7', 
    'Reporter intensity corrected 5 K3971(126).K3972(127).K39711(128).K3975(130).K39710(131).R7', 
    'Reporter intensity corrected 6 K3971(126).K3972(127).K39711(128).K3975(130).K39710(131).R7', 
    
    'Reporter intensity corrected 1 K3971(126).K3972(127).K3973(128).K3974(129).BLANK(130).K3976(131).R5', 
    'Reporter intensity corrected 2 K3971(126).K3972(127).K3973(128).K3974(129).K3976(131).R5', 
    'Reporter intensity corrected 3 K3971(126).K3972(127).K3973(128).K3974(129).K3976(131).R5', 
    'Reporter intensity corrected 4 K3971(126).K3972(127).K3973(128).K3974(129).K3976(131).R5', 
    'Reporter intensity corrected 5 K3971(126).K3972(127).K3973(128).K3974(129).K3976(131).R5', 
    'Reporter intensity corrected 6 K3971(126).K3972(127).K3973(128).K3974(129).K3976(131).R5', 
    
    'Reporter intensity corrected 1 K3971(126).K3972(127).K3977(128).K3978(129).K3979(130).R6', 
    'Reporter intensity corrected 2 K3971(126).K3972(127).K3977(128).K3978(129).K3979(130).R6', 
    'Reporter intensity corrected 3 K3971(126).K3972(127).K3977(128).K3978(129).K3979(130).R6', 
    'Reporter intensity corrected 4 K3971(126).K3972(127).K3977(128).K3978(129).K3979(130).R6', 
    'Reporter intensity corrected 5 K3971(126).K3972(127).K3977(128).K3978(129).K3979(130).R6', 
    'Reporter intensity corrected 6 K3971(126).K3972(127).K3977(128).K3978(129).K3979(130).R6', 
    
    'Reporter intensity corrected 1 K3981(126).K3982(127).K3983(128).K3984(129).R8', 
    'Reporter intensity corrected 2 K3981(126).K3982(127).K3983(128).K3984(129).R8', 
    'Reporter intensity corrected 3 K3981(126).K3982(127).K3983(128).K3984(129).R8', 
    'Reporter intensity corrected 4 K3981(126).K3982(127).K3983(128).K3984(129).R8', 
    'Reporter intensity corrected 5 K3981(126).K3982(127).K3983(128).K3984(129).R8', 
    'Reporter intensity corrected 6 K3981(126).K3982(127).K3983(128).K3984(129).R8', 
    
    'Reporter intensity corrected 1 K3991(126).K3992(127).K3993(128).K3994(129).R9', 
    'Reporter intensity corrected 2 K3991(126).K3992(127).K3993(128).K3994(129).R9', 
    'Reporter intensity corrected 3 K3991(126).K3992(127).K3993(128).K3994(129).R9', 
    'Reporter intensity corrected 4 K3991(126).K3992(127).K3993(128).K3994(129).R9', 
    'Reporter intensity corrected 5 K3991(126).K3992(127).K3993(128).K3994(129).R9', 
    'Reporter intensity corrected 6 K3991(126).K3992(127).K3993(128).K3994(129).R9'
)
demultiplex(dequantify(x))
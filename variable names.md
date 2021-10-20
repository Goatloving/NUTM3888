

# Dataset

## Dat2



# Y

## Y1

- all responses (xxx_scores) are binary (0,1)

## Y2


# X

## X1: 

- all predictors are encoded as categorical variables: Sys_score, Dis_score,Tri_score, Chol_score, LDL_score, Glu_score, HDL_score, Waist_score, MBC_score, VBC_score, ApoB_score, HbA1c_score
- no BMI
- no SMOKE
- exerices encoded "wrongly" (0-> longer time, 1->shorter time)

## X2:

- keep numerical values
- categorical values encoeded as X1


# Models

## Model 1:

Y1$cvd_score~X1

## Model 2:


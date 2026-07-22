# Model Tuning workflow

## Fixed setting

* NPE
* Clonal tree
* Training data: 10000
* Testing data: 1000

## Performance measures

* Validation negative log likelihood (main)
* SBC p-values
* Distance between posterior mean and true values

## Hyperparameters and detailed methods

### Input Dim

Always exclude LD and G4 measures under 100 bp.

* Input dim: full, no exact LD and G4, no prop LD and G4

### Summary NN

* Hidden layers: 1, 2, 3, 4
* Hidden units: 32, 48, 64, 128, 256
* Output dim: $2+1, 2\times2+1, 2\times4+1, 2\times8+1$

### Normalizing flow

* Density estimator: NSF, MAF
* Number of transforms: 3, 5, 8
* Hidden features: 30, 50, 80, 120
* Spline bins (NSF): 5, 10, 20

### NPE training

* Learning rate: 0.0001, 0.0005, 0.001
* Weight decay: $0, 10^{-6}, 10^{-5}, 10^{-4}$
* Stop after epochs: 10, 20, 30
* Training batch size: 100, 200, 300
* Clip max norm: None, 1, 5, 8

## Stage 0

Basic model with no summary NN:
* Density estimator: NSF
* Number of transforms: 5
* Hidden features: 50
* Spline bins (NSF): 10
* Learning rate: 0.0005
* Weight decay: $0$
* Stop after epochs: 20
* Training batch size: 200
* Clip max norm: 5

**DONE**

## Stage 1

Baseline model config, use the default setting for normalizing flow and training:
* Input dim: full
* Hidden layers: 2
* Hidden units: 48
* Output dim: $2\times2+1$
* Density estimator: NSF
* Number of transforms: 5
* Hidden features: 50
* Spline bins (NSF): 10
* Learning rate: 0.0005
* Weight decay: $0$
* Stop after epochs: 20
* Training batch size: 200
* Clip max norm: 5

**DONE**


## Stage 2

Input dim

Find the best one, if all simular, choose full dim.

**DONE - choose full dim**

## Stage 3

Output dim

Find the best two, pass to the next stage.

**DONE - choose $2+1, 2\times2+1$**

## Stage 4

Hidden layers and units

Broad search: 1 run; top 8 find top 3: 5 runs. Find the best three, pass to the next stage.

**DONE - best three:**
1. Output: 2+1, layers: 1, units: 128
2. Output: 2+1, layers: 4, units: 256
3. Output: 2+1, layers: 4, units: 48

## Stage 5

NSF or MAF

Find the best one.

**DONE - choose NSF**

## Stage 6

Config for normalizing flow: number of transforms, hidden features, spline bins

Broad search: 1 run; top 8 find top 3: 5 runs

## Stage 7

NPE training setting: learning rate, weight decay, stop after epochs, training batch size, clip max norm

Broad search: 1 run; top 8 find top 3: 5 runs

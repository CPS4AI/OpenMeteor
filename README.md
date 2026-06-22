# OpenMeteor

OpenMeteor is a semi-honest 3-party neural-network inference framework based on the Meteor protocol: masked fixed-point arithmetic, replicated secret sharing, setup-generated correlated randomness, and an explicit split between preprocessing and online execution costs.

This repository is a research and benchmarking codebase. It is not a production cryptographic library, but it now runs local secure inference end to end on the bundled preloaded MNIST networks.

## What Is Inside

- Secure 3PC inference primitives for FC, CNN, BN, ReLU, Maxpool, MatMul, DotProduct, PrivateCompare, Bit2A, PreMSB support material, and faithful fixed-point truncation.
- Setup/preprocessing generation for the correlated randomness consumed by the active online protocol path.
- Preloaded MNIST network shares for `SecureML`, `Sarda`, `MiniONN`, and `LeNet`.
- Per-run cost reports split into:
  - preprocessing computation time
  - online computation time
  - preprocessing communication
  - online communication
- Dual build path:
  - native arm64 portable AES backend
  - x86_64 AES-NI backend, including Rosetta on Apple Silicon

## Repository Layout

```text
files/          Party keys, IP files, and bundled preloaded network shares
files/preload/  Preloaded MNIST inputs, weights, biases, references, and predictions
lib_eigen/      Eigen matrix library
src/            Protocols, layers, communication, setup, and runners
util/           AES/hash/GF helper code
scripts/        Model/share generation notebooks and scripts
docs/           Implementation notes and protocol audit trail
```

## Build

OpenMeteor uses `make` and a C++11 compiler.

On Apple Silicon or other arm64 hosts, build the native portable backend:

```bash
make -B all TARGET_ARCH=native -j2
```

On x86_64, or on Apple Silicon through Rosetta, build the AES-NI backend:

```bash
make -B all TARGET_ARCH=x86_64 -j2
```

Inspect the selected backend:

```bash
make backend TARGET_ARCH=native
make backend TARGET_ARCH=x86_64
```

Linux x86_64 builds may need OpenSSL development headers:

```bash
sudo apt-get install g++ make libssl-dev
```

## Run

Run a local 3-party unit test:

```bash
make command \
  RUN_MODE=unit \
  NETWORK=SecureML \
  DATASET=MNIST \
  UNIT_TEST=MeteorTruncation \
  PORT_BASE=38100
```

Run preloaded secure inference:

```bash
make command \
  RUN_MODE=preloaded \
  NETWORK=SecureML \
  DATASET=MNIST \
  UNIT_TEST=MeteorRELU \
  PORT_BASE=38200
```

Useful variables:

```text
TARGET_ARCH  native, arm64, or x86_64
RUN_MODE     unit, inference, preloaded, or train
NETWORK      SecureML, Sarda, MiniONN, LeNet, AlexNet, or VGG16
DATASET      MNIST, CIFAR10, or ImageNet
UNIT_TEST    MeteorTruncation, MeteorRELU, MeteorPC, MeteorMaxpool, Conv, BN, ...
PORT_BASE    first localhost port used by the 3 local parties
```

## Cost Reports

Every measured run prints local and aggregated costs. Party 0 aggregates communication across all three parties and reports preprocessing and online buckets separately:

```text
Aggregated cost breakdown for SecureML preloaded:
Preprocessing computation: 0.0682206 sec max wall, 0.058805 sec total CPU
Online computation: 0.0909841 sec max wall, 0.226087 sec total CPU

Total communication: 3.33399MB (sent) and 3.33399MB (recv)
Preprocessing communication: 3.23217MB (sent) and 3.23217MB (recv), 564 sends and 564 recvs
Online communication: 0.101826MB (sent) and 0.101826MB (recv), 105 sends and 105 recvs
```

Preprocessing includes setup-generation calls such as random masks, product correlations, Bit2A material, PreMSB material, PrivateCompare material, and faithful truncation masks. Online includes the protocol work that consumes those correlations and opens masked values.

## Verified Local Targets

Recent local checks include:

```bash
make -B all TARGET_ARCH=native -j2
make -B all TARGET_ARCH=x86_64 -j2
make command RUN_MODE=unit NETWORK=SecureML DATASET=MNIST UNIT_TEST=MeteorTruncation PORT_BASE=38100 TARGET_ARCH=native
make command RUN_MODE=unit NETWORK=SecureML DATASET=MNIST UNIT_TEST=MeteorTruncation PORT_BASE=38200 TARGET_ARCH=x86_64
make command RUN_MODE=preloaded NETWORK=SecureML DATASET=MNIST UNIT_TEST=MeteorRELU PORT_BASE=37900 TARGET_ARCH=x86_64
make command RUN_MODE=preloaded NETWORK=SecureML DATASET=MNIST UNIT_TEST=MeteorRELU PORT_BASE=38000 TARGET_ARCH=native
```

The bundled preloaded networks have also been checked against their reference final-output argmax values:

```text
SecureML -> 6
Sarda    -> 1
MiniONN  -> 1
LeNet    -> 8
```

## Scope

OpenMeteor targets semi-honest 3PC only and is protoype-of-concept, which should not be used in product. The arm64 backend is designed for portability and development; the x86_64 AES-NI backend remains the performance-oriented path.

## Citation

```text
@inproceedings{10.1145/3543507.3583272,
author = {Dong, Ye and Xiaojun, Chen and Jing, Weizhan and Kaiyun, Li and Wang, Weiping},
title = {Meteor: Improved Secure 3-Party Neural Network Inference with Reducing Online Communication Costs},
year = {2023},
isbn = {9781450394161},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3543507.3583272},
doi = {10.1145/3543507.3583272},
booktitle = {Proceedings of the ACM Web Conference 2023},
pages = {2087-2098},
location = {Austin, TX, USA},
series = {WWW '23}
}
```

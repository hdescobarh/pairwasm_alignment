# pairwasm_alignment

![Rust](https://img.shields.io/badge/-Rust-B7410E?logo=rust&logoColor=28282B&labelColor=white)
![WebAssembly](https://img.shields.io/badge/-WebAssembly-654FF0?logo=webassembly&logoColor=654FF0&labelColor=white)
![npm](https://img.shields.io/badge/-npm-CC3534?logo=npm&labelColor=gray)
![Experimental](https://img.shields.io/badge/stability-experimental-orange)
![License](https://img.shields.io/badge/license-MIT-blue)

Leveraging the performance, portability, and interoperability of WebAssembly (Wasm) brings the advantages of cloud computing to bioinformatics and computational biology, enabling seamless collaboration among researchers and enhancing reproducibility. Wasm runs efficiently across various operating systems and hardware architectures, delivering near-native performance without the overhead of traditional virtual machines. Its low system requirements also enable rapid scaling, making cloud computing more accessible and cost-effective.

This project introduces an experimental Wasm-based module for pairwise sequence alignment, developed in Rust 🦀. It implements the classic Needleman-Wunsch and Smith-Waterman algorithms, commonly used to align DNA, RNA, and protein sequences. By combining Wasm’s cross-platform efficiency with Rust’s performance and safety, this module offers a portable, high-speed solution for sequence alignment, ideal for both cloud and web applications.

## Crate documentation and demonstration

The WASM module use is straightforward since it only exposes a single function; read the next section for details. Also, you can check the Rust crate [documentation here](https://hdescobarh.github.io/pairwasm_alignment/pairwasm_alignment/).

I also deployed an online demonstration that runs __locally__ and on __your browser__. Be aware that Smith-Waterman and Needleman-Wunsch are dynamic programming algorithms with [quadratic time complexity](https://en.wikipedia.org/wiki/Big_O_notation). __Do not__ use the demo for long sequences.

</br>

[__Check the demo__ </br>![pairwasm_alignment](image.png)](https://hansescobar.com/en/demo/pairwasm)

## Releases and how to use them

The module is offered in two flavors:

- [npm Module](https://www.npmjs.com/package/pairwasm_alignment)

The wasm module itself is natively an ES module [^1]. It needs a Bundler. Only Webpack offers total compatibility and other bundlers will require additional configurations.

- [Web](https://github.com/hdescobarh/pairwasm_alignment/releases)

It can natively be included on a web page, and doesn't require any further postprocessing [^1].

For example,

```typescript
import init, { do_protein_alignment } from "./wasm_module/pairwasm_alignment.js";

export async function run(
  string_1: string,
  string_2: string,
  open_cost: number,
  extend_cost: number,
  substitution_matrix: number,
  algorithm: number) {
  await init();
  return do_protein_alignment(
    string_1,
    string_2,
    open_cost,
    extend_cost,
    substitution_matrix,
    algorithm
  );
}
```

Then the *run* function can be imported normally in any script.

## Known issues

- Using local alignment can give suboptimal alignments.
</br></br>

--------

[^1]: [Deploying Rust and WebAssembly](https://rustwasm.github.io/docs/wasm-bindgen/reference/deployment.html#deploying-rust-and-webassembly).

## Bibliography

- Chao, J., Tang, F., & Xu, L. (2022). Developments in Algorithms for Sequence Alignment: A Review. Biomolecules, 12(4), 546. <https://doi.org/10.3390/biom12040546>
- Masek, W. J., & Paterson, M. S. (1980). A faster algorithm computing string edit distances. Journal of Computer and System Sciences, 20(1), 18–31. <https://doi.org/10.1016/0022-0000(80)90002-1>
- Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. Journal of Molecular Biology, 48(3), 443–453. <https://doi.org/10.1016/0022-2836(70)90057-4>
- Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences. Journal of Molecular Biology, 147(1), 195–197. <https://doi.org/10.1016/0022-2836(81)90087-5>
- Smith, T. F., Waterman, M. S., & Fitch, W. M. (1981). Comparative biosequence metrics. Journal of Molecular Evolution, 18(1), 38–46. <https://doi.org/10.1007/BF01733210>

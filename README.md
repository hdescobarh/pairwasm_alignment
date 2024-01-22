# pairwasm_alignment

![Rust](https://img.shields.io/badge/-Rust-B7410E?logo=rust&logoColor=28282B&labelColor=white)
![WebAssembly](https://img.shields.io/badge/-WebAssembly-654FF0?logo=webassembly&logoColor=654FF0&labelColor=white)
![npm](https://img.shields.io/badge/-npm-CC3534?logo=npm&labelColor=white)
![Experimental](https://img.shields.io/badge/stability-experimental-orange)
![License](https://img.shields.io/badge/license-MIT-blue)

## Description

This is an experiment implementing a WebAssembly pairwise sequence alignment module in Rust 🦀.

> [!NOTE]
> [Rust crate documentation](https://hdescobarh.github.io/pairwasm_alignment/pairwasm_alignment/)

## Releases

- [npm Module](https://www.npmjs.com/package/pairwasm_alignment)

> The wasm module itself is natively an ES module. It needs a Bundler [^1].

- [Web]

> It can natively be included on a web page, and doesn't require any further postprocessing [^1].

## Known issues

- Using affine gap model can give suboptimal alignments.

[^1]: [Deploying Rust and WebAssembly](https://rustwasm.github.io/docs/wasm-bindgen/reference/deployment.html#deploying-rust-and-webassembly).

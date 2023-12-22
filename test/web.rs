//! Test suite for the Web and headless browsers.
#![cfg(target_arch = "wasm32")]

// https://rustwasm.github.io/wasm-bindgen/wasm-bindgen-test/usage.html
// the tests must be in the root of the crate, or within a pub mod
use wasm_bindgen_test::*;

// This macro macro is used to indicate that the test is intended to be
// executed in a web browser (and not in other environment; Node.js for example)
wasm_bindgen_test_configure!(run_in_browser);

#[wasm_bindgen_test]
fn pass() {
    let x = 5 - 3;
    assert_eq!(x, 2);
}

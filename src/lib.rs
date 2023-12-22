pub mod bioseq;
pub mod matrix;
pub mod scoring;
pub mod utils;

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
}

#[wasm_bindgen]
pub fn greet() {
    alert("Hello, {{project-name}}!");
}

pub fn add() -> usize {
    1 + 1
}
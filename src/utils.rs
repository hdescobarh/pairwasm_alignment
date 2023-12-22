// https://rustwasm.github.io/docs/wasm-pack/tutorials/npm-browser-packages/template-deep-dive/src-utils-rs.html
pub fn set_panic_hook() {
    // When the `console_error_panic_hook` feature is enabled, we can call the
    // `set_panic_hook` function at least once during initialization, and then
    // we will get better error messages if our code ever panics.
    //
    // For more details see
    // https://github.com/rustwasm/console_error_panic_hook#readme
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

/// Flag to denote the minimal units of a sequence.
pub trait AlignmentUnit
where
    Self: Sized + Ord + Eq + Copy,
{
    /// Represents each Aac duple with a unique integer identifier.
    ///
    /// It is a map (amino acid code, amino acid code) ↦ integer.\n
    /// The current implementation is the cantor pairing function,
    /// which is bijective and strictly monotonic.
    fn duple_pairing(code_1: Self, code_2: Self) -> u16;

    /// Cantor pairing function.
    /// A bijective and strictly monotonic map (𝒂, 𝒃) ↦ 𝒄, such that 𝒂,𝒃,𝒄 ∊ ℕ₀
    fn cantor_pairing(k1: u16, k2: u16) -> u16 {
        let k12 = k1 + k2;
        (k12 * (k12 + 1)).div_euclid(2) + k2
    }
}

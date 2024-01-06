//! Deals with the output format

use crate::aligner::utils::AlignmentSequence;
use crate::bioseq::Aac;
use crate::utils::AlignmentUnit;
use std::cmp::PartialEq;
use std::convert::From;
use std::fmt::Display;

const GAP_STR: char = '_';
const MATCH_STR: char = '|';
const MISMATCH_STR: char = ':';
const SPACE_STR: char = '\u{0020}';

impl From<&Aac> for char {
    fn from(val: &Aac) -> Self {
        match val {
            Aac::A => 'A',
            Aac::C => 'C',
            Aac::D => 'D',
            Aac::E => 'E',
            Aac::F => 'F',
            Aac::G => 'G',
            Aac::H => 'H',
            Aac::I => 'I',
            Aac::K => 'K',
            Aac::L => 'L',
            Aac::M => 'M',
            Aac::N => 'N',
            Aac::P => 'P',
            Aac::Q => 'Q',
            Aac::R => 'R',
            Aac::S => 'S',
            Aac::T => 'T',
            Aac::V => 'V',
            Aac::W => 'W',
            Aac::Y => 'Y',
        }
    }
}

impl<A> Display for AlignmentSequence<A>
where
    A: AlignmentUnit + PartialEq,
    char: for<'a> From<&'a A>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut line1: String = String::with_capacity(self.read().len());
        let mut line2: String = String::with_capacity(self.read().len());
        let mut line3: String = String::with_capacity(self.read().len());
        let mut line_break_counter = 0;
        for [first, second] in self.read() {
            if line_break_counter == 50 {
                line1.push('\n');
                line2.push('\n');
                line3.push('\n');
                line_break_counter = 0;
            }
            match (first, second) {
                (None, None) => {
                    line1.push(GAP_STR);
                    line2.push(SPACE_STR);
                    line3.push(GAP_STR);
                }
                (None, Some(aa)) => {
                    line1.push(GAP_STR);
                    line2.push(SPACE_STR);
                    line3.push(aa.into());
                }
                (Some(aa), None) => {
                    line1.push(aa.into());
                    line2.push(SPACE_STR);
                    line3.push(GAP_STR);
                }
                (Some(aa1), Some(aa2)) => {
                    if aa1 == aa2 {
                        line2.push(MATCH_STR);
                    } else {
                        line2.push(MISMATCH_STR);
                    }
                    line1.push(aa1.into());

                    line3.push(aa2.into());
                }
            };
        }

        f.write_fmt(format_args!("{}\n{}\n{}", line1, line2, line3))
    }
}

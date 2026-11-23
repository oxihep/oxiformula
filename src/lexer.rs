//! Tokenizer for ROOT `TFormula` expressions. Every token keeps the byte range
//! of the source it came from, so errors can point at it.

use std::ops::Range;

use crate::ParseError;

/// A lexical token.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Tok {
    /// A numeric literal (`1.5`, `2e-3`, `.5`).
    Num(f64),
    /// An identifier, `::` included: a variable, a constant, a shortcut
    /// (`gaus`, `pol2`) or a function (`sin`, `TMath::Gaus`).
    Ident(String),
    /// The trimmed contents of `[…]`: a parameter (`0`, `p1`, `mean`), the
    /// index of `x[1]`, or a shortcut's variables (`gaus[x,y]`).
    Bracket(String),
    Plus,
    Minus,
    Star,
    Slash,
    /// Exponentiation, from `^` or `**`.
    Caret,
    LParen,
    RParen,
    Comma,
    Lt,
    Gt,
    Le,
    Ge,
    EqEq,
    Ne,
    And,
    Or,
    Not,
    Question,
    Colon,
}

#[derive(Debug, Clone)]
pub(crate) struct Token {
    pub(crate) tok: Tok,
    pub(crate) span: Range<usize>,
}

/// Turn `src` into a token stream, or a [`ParseError`] on a character that
/// cannot start a token.
pub(crate) fn lex(src: &str) -> Result<Vec<Token>, ParseError> {
    let bytes = src.as_bytes();
    let mut i = 0;
    let mut out = Vec::new();
    let next = |i: usize| bytes.get(i + 1).copied();
    while i < bytes.len() {
        let start = i;
        let (tok, len) = match bytes[i] {
            b' ' | b'\t' | b'\n' | b'\r' => {
                i += 1;
                continue;
            }
            b'+' => (Tok::Plus, 1),
            b'-' => (Tok::Minus, 1),
            b'*' if next(i) == Some(b'*') => (Tok::Caret, 2),
            b'*' => (Tok::Star, 1),
            b'/' => (Tok::Slash, 1),
            b'^' => (Tok::Caret, 1),
            b'(' => (Tok::LParen, 1),
            b')' => (Tok::RParen, 1),
            b',' => (Tok::Comma, 1),
            b'<' if next(i) == Some(b'=') => (Tok::Le, 2),
            b'<' => (Tok::Lt, 1),
            b'>' if next(i) == Some(b'=') => (Tok::Ge, 2),
            b'>' => (Tok::Gt, 1),
            b'=' if next(i) == Some(b'=') => (Tok::EqEq, 2),
            b'!' if next(i) == Some(b'=') => (Tok::Ne, 2),
            b'!' => (Tok::Not, 1),
            b'&' if next(i) == Some(b'&') => (Tok::And, 2),
            b'|' if next(i) == Some(b'|') => (Tok::Or, 2),
            b'?' => (Tok::Question, 1),
            b':' if next(i) == Some(b':') => {
                return Err(ParseError::new("`::` must follow a name", i..i + 2)
                    .with_hint("write a namespaced function as `TMath::Gaus(…)`"));
            }
            b':' => (Tok::Colon, 1),
            b'[' => {
                let close = src[i + 1..]
                    .find(']')
                    .ok_or_else(|| ParseError::new("`[` is never closed", i..i + 1))?;
                let inner = src[i + 1..i + 1 + close].trim();
                if inner.is_empty() {
                    return Err(ParseError::new("empty brackets `[]`", i..i + close + 2)
                        .with_hint("a parameter is `[0]`, `[1]`, … or a name such as `[mean]`"));
                }
                (Tok::Bracket(inner.to_owned()), close + 2)
            }
            b'0'..=b'9' | b'.' => number(src, i)?,
            b'A'..=b'Z' | b'a'..=b'z' | b'_' => {
                let mut j = i;
                while j < bytes.len() {
                    match bytes[j] {
                        b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'_' => j += 1,
                        // `TMath::Sin` is one name.
                        b':' if bytes.get(j + 1) == Some(&b':')
                            && bytes
                                .get(j + 2)
                                .is_some_and(|c| c.is_ascii_alphabetic() || *c == b'_') =>
                        {
                            j += 2;
                        }
                        _ => break,
                    }
                }
                (Tok::Ident(src[i..j].to_owned()), j - i)
            }
            _ => return Err(illegal(src, i)),
        };
        out.push(Token {
            tok,
            span: start..start + len,
        });
        i = start + len;
    }
    Ok(out)
}

/// A numeric literal starting at `start`: digits, one `.`, and an exponent.
fn number(src: &str, start: usize) -> Result<(Tok, usize), ParseError> {
    let bytes = src.as_bytes();
    let mut i = start;
    let (mut dot, mut exp) = (false, false);
    while i < bytes.len() {
        match bytes[i] {
            b'0'..=b'9' => i += 1,
            b'.' if !dot && !exp => {
                dot = true;
                i += 1;
            }
            b'e' | b'E' if !exp => {
                exp = true;
                i += 1;
                if matches!(bytes.get(i), Some(b'+' | b'-')) {
                    i += 1;
                }
            }
            _ => break,
        }
    }
    let text = &src[start..i];
    let value = text
        .parse::<f64>()
        .map_err(|_| ParseError::new(format!("`{text}` is not a number"), start..i))?;
    if !value.is_finite() {
        return Err(
            ParseError::new(format!("`{text}` is too large for a double"), start..i)
                .with_hint("write infinity as `TMath::Infinity()`"),
        );
    }
    Ok((Tok::Num(value), i - start))
}

/// The error for a character that starts no token, with a hint for the
/// spellings other languages use.
fn illegal(src: &str, i: usize) -> ParseError {
    let c = src[i..].chars().next().unwrap_or('?');
    let span = i..i + c.len_utf8();
    let e = ParseError::new(format!("unexpected character `{c}`"), span);
    match c {
        '&' => e.with_hint("ROOT writes logical AND as `&&`"),
        '|' => e.with_hint("ROOT writes logical OR as `||`"),
        '~' => e.with_hint("ROOT writes logical NOT as `!`"),
        '=' => e.with_hint("compare with `==`"),
        '%' => e.with_hint("TFormula compiles to C++, where `%` takes integers: use `fmod(a, b)`"),
        _ => e,
    }
}

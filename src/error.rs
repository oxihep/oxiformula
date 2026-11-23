//! The error a formula that does not parse gives.

use std::fmt;
use std::ops::Range;

/// A formula that does not parse: what is wrong, where (a byte range of the
/// source), and, for the mistakes people make most often, how to fix it.
///
/// Its `Display` shows the formula with the offending part underlined:
///
/// ```text
/// unknown function `bogus`
///   [0]*bogus(x)
///       ^^^^^
///   help: the functions are those of ROOT's TFormula (sin, exp, TMath::Gaus, …)
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    message: String,
    span: Range<usize>,
    hint: Option<String>,
    source: String,
}

impl ParseError {
    pub(crate) fn new(message: impl Into<String>, span: Range<usize>) -> ParseError {
        ParseError {
            message: message.into(),
            span,
            hint: None,
            source: String::new(),
        }
    }

    pub(crate) fn with_hint(mut self, hint: impl Into<String>) -> ParseError {
        self.hint = Some(hint.into());
        self
    }

    /// Attach the formula the span points into, for `Display`.
    pub(crate) fn in_source(mut self, source: &str) -> ParseError {
        source.clone_into(&mut self.source);
        let end = source.len();
        self.span = self.span.start.min(end)..self.span.end.min(end);
        self
    }

    /// What is wrong, in one line.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    /// The byte offset into the formula where the problem starts.
    #[must_use]
    pub fn position(&self) -> usize {
        self.span.start
    }

    /// The byte range of the formula the problem covers (empty at the end of
    /// a formula that stops too early).
    #[must_use]
    pub fn span(&self) -> Range<usize> {
        self.span.clone()
    }

    /// How to fix it, when the mistake is a common one.
    #[must_use]
    pub fn hint(&self) -> Option<&str> {
        self.hint.as_deref()
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)?;
        if !self.source.is_empty() && !self.source.contains('\n') {
            let before = self.source[..self.span.start].chars().count();
            let width = self.source[self.span.clone()].chars().count().max(1);
            let line: String = self
                .source
                .chars()
                .map(|c| if c.is_control() { ' ' } else { c })
                .collect();
            write!(
                f,
                "\n  {line}\n  {}{}",
                " ".repeat(before),
                "^".repeat(width)
            )?;
        } else if self.source.contains('\n') {
            write!(f, " (at byte {})", self.span.start)?;
        }
        if let Some(hint) = &self.hint {
            write!(f, "\n  help: {hint}")?;
        }
        Ok(())
    }
}

impl std::error::Error for ParseError {}

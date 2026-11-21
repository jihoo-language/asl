use std::fmt;

#[derive(Clone, Copy)]
pub(crate) struct Span {
    pub(crate) line: usize,
    pub(crate) column: usize,
}

pub(crate) struct Error {
    pub(crate) span: Option<Span>,
    pub(crate) message: String,
}

impl Error {
    pub(crate) fn plain(message: impl Into<String>) -> Self {
        Self {
            span: None,
            message: message.into(),
        }
    }

    pub(crate) fn new(span: Span, message: impl Into<String>) -> Self {
        Self::plain(message).at(span)
    }

    pub(crate) fn at(mut self, span: Span) -> Self {
        self.span.get_or_insert(span);
        self
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.span {
            Some(span) => write!(f, ":{}:{}: {}", span.line, span.column, self.message),
            None => write!(f, ": {}", self.message),
        }
    }
}

pub(crate) type Result<T> = std::result::Result<T, Error>;

use crate::error::{Error, Result, Span};
use crate::map_bracket::Operation as MapBracketOp;
use crate::map_system::Operation as MapSystemOp;
use crate::map_transfer::Operation as MapTransferOp;
use crate::odd_cobar::FactorLiteral;
use crate::odd_massey::Operation as OddSystemOp;
use crate::transfer::Operation as TransferOp;

const MAX_NESTING: usize = 64;
const MAX_IDENTIFIER_BYTES: usize = 128;

#[derive(Clone)]
pub(crate) struct Label {
    pub(crate) text: String,
    pub(crate) span: Span,
}

pub(crate) struct ParsedType {
    pub(crate) span: Span,
    pub(crate) kind: TypeForm,
    depth: usize,
}

pub(crate) enum TypeForm {
    Name(Label),
    Star(Box<ParsedType>, Label),
    Cohomology(Box<ParsedType>, usize),
    Operations(usize),
    Cobar(usize, usize),
    PCobar(usize, usize, usize),
    PCobarAn(usize, usize, usize, usize),
    PExt(usize, usize, usize),
    PExtAn(usize, usize, usize, usize),
    PMay(usize, usize, usize, usize, usize),
    PMayAn(usize, usize, usize, usize, usize, usize),
    CobarAn(usize, usize, usize),
    Ext(usize, usize),
    ExtAn(usize, usize, usize),
    XExt(Box<ParsedType>, usize, usize),
    MapCobar(Box<ParsedType>, Box<ParsedType>, usize, i64),
    MapExt(Box<ParsedType>, Box<ParsedType>, usize, i64),
    MapMay(MapMaySpec),
    QKoszul(QSpec),
    QExt(QSpec),
    VKoszul(QSpec<i64>),
    VExt(QSpec<i64>),
    VBock(VBockSpec),
    May(usize, usize, usize, usize),
    MayAn(usize, usize, usize, usize, usize),
    CeAn(usize, usize, usize, usize, usize),
    Suspension(Box<ParsedType>),
    Stable(Box<ParsedType>, usize, usize),
    Bar(usize),
    CyclicH2(usize),
    Product(Box<ParsedType>, Box<ParsedType>),
    Arrow(Box<ParsedType>, Box<ParsedType>),
}

pub(crate) struct QSpec<S = usize> {
    pub(crate) source: Box<ParsedType>,
    pub(crate) target: Box<ParsedType>,
    pub(crate) selected: Vec<usize>,
    pub(crate) s: S,
    pub(crate) t: i64,
}

pub(crate) struct VBockSpec {
    pub(crate) localized: QSpec<i64>,
    pub(crate) r: usize,
    pub(crate) f: usize,
}

pub(crate) struct MapMaySpec {
    pub(crate) source: Box<ParsedType>,
    pub(crate) target: Box<ParsedType>,
    pub(crate) r: usize,
    pub(crate) s: usize,
    pub(crate) t: i64,
    pub(crate) f: usize,
}

pub(crate) struct ParsedExpr {
    pub(crate) span: Span,
    pub(crate) kind: ExprForm,
    depth: usize,
}

pub(crate) type Entries = Vec<(Label, ParsedExpr)>;

#[derive(Clone, Copy)]
pub(crate) enum BinaryOp {
    Add,
    Cup,
    Pull,
    Mul,
    Act,
    Concat,
    Yoneda,
    XAction,
    MapExtCompose,
    XMap,
    MayMul,
    CeMul,
    MilnorPair,
}

impl BinaryOp {
    fn from_name(name: &str) -> Option<Self> {
        match name {
            "add" => Some(Self::Add),
            "cup" => Some(Self::Cup),
            "pull" => Some(Self::Pull),
            "mul" => Some(Self::Mul),
            "act" => Some(Self::Act),
            "concat" => Some(Self::Concat),
            "yoneda" => Some(Self::Yoneda),
            "xaction" => Some(Self::XAction),
            "ecomp" => Some(Self::MapExtCompose),
            "xmap" => Some(Self::XMap),
            "mayMul" => Some(Self::MayMul),
            "ceMul" => Some(Self::CeMul),
            "milnorPair" => Some(Self::MilnorPair),
            _ => None,
        }
    }
}

pub(crate) enum ClassForm {
    Zero,
    Generator(usize, Span),
    Cocycle(Vec<Vec<Label>>),
}

pub(crate) enum SuspensionForm {
    North,
    South,
    Equator,
}

pub(crate) enum CyclicForm {
    Support(Vec<(usize, Span)>),
    Binomial(usize, Span),
    Carry,
}

pub(crate) struct Witness {
    pub(crate) from: usize,
    pub(crate) to: usize,
    pub(crate) span: Span,
    pub(crate) value: ParsedExpr,
}

pub(crate) fn validate_witness_interval(
    arity: usize,
    from: usize,
    to: usize,
    span: Span,
) -> Result<()> {
    if from >= to || to > arity {
        return Err(Error::new(span, format!("invalid interval [{from},{to}]")));
    }
    if to - from == 1 {
        return Err(Error::new(
            span,
            format!("adjacent interval [{from},{to}] already has an input"),
        ));
    }
    if to - from == arity {
        return Err(Error::new(
            span,
            format!("top interval [{from},{to}] cannot have a witness"),
        ));
    }
    Ok(())
}

pub(crate) struct QTerm {
    pub(crate) source: MapGenerator,
    pub(crate) target: MapGenerator,
    pub(crate) selected: Vec<usize>,
    pub(crate) powers: Vec<usize>,
}

pub(crate) struct MapGenerator {
    pub(crate) base: ParsedType,
    pub(crate) degree: usize,
    pub(crate) index: usize,
    pub(crate) span: Span,
}

pub(crate) enum ExprForm {
    Name(Label),
    Point(ParsedType, Label),
    Class(ParsedType, usize, ClassForm),
    Zero(ParsedType),
    Square(usize),
    Milnor(Vec<usize>),
    Coact(usize, usize, Span, Box<ParsedExpr>),
    CobarWord(Vec<(Vec<usize>, Span)>),
    PCobarWord(usize, Vec<FactorLiteral>),
    PExtGenerator(usize, usize, usize, usize, Span),
    PExtAnGenerator(usize, usize, usize, usize, usize, Span),
    PMayGenerator(usize, usize, usize, usize, usize, usize, Span),
    PMayAnGenerator(usize, usize, usize, usize, usize, usize, usize, Span),
    PScale(usize, Box<ParsedExpr>),
    ExtGenerator(usize, usize, usize, Span),
    ExtAnGenerator(usize, usize, usize, usize, Span),
    RestrictAlgebra(usize, Box<ParsedExpr>),
    CobarCup(usize, Box<ParsedExpr>, Box<ParsedExpr>),
    ExtSquare(usize, Box<ParsedExpr>),
    Nakamura(usize, Box<ParsedExpr>),
    MayLiftTo(usize, usize, Box<ParsedExpr>),
    HigherRestriction(usize, Vec<ParsedExpr>),
    XExtGenerator(ParsedType, usize, usize, usize, Span),
    MapWord(MapGenerator, MapGenerator, Box<ParsedExpr>),
    MapExtGenerator(ParsedType, ParsedType, usize, i64, usize, Span),
    MapExtIdentity(ParsedType),
    MapMayGenerator(MapMaySpec, usize, Span),
    QExtGenerator(QSpec, usize, Span),
    QTerm(QTerm),
    QIdentity(ParsedType, Vec<usize>),
    VExtGenerator(QSpec<i64>, usize, Span),
    VBockGenerator(VBockSpec, usize, Span),
    VBockClass(usize, usize, Box<ParsedExpr>),
    VShift(i64, Box<ParsedExpr>),
    QRestrict(Vec<usize>, Box<ParsedExpr>),
    MayClass(usize, usize, Box<ParsedExpr>),
    MayGenerator(usize, usize, usize, usize, usize, Span),
    MayAnGenerator(usize, usize, usize, usize, usize, usize, Span),
    CeGenerator(usize, usize, usize, usize, usize, usize, Span),
    CeClass(usize, usize, Box<ParsedExpr>),
    Suspension(ParsedType, SuspensionForm),
    Shift(usize, Box<ParsedExpr>),
    Cyclic(usize, Span, CyclicForm),
    Massey(Vec<ParsedExpr>, Vec<Witness>),
    PowerDef(Vec<ParsedExpr>),
    Bracket(Vec<ParsedExpr>),
    Transfer(TransferOp, Vec<ParsedExpr>),
    MapBracket(MapBracketOp, Vec<ParsedExpr>),
    MapSystem(MapSystemOp, Vec<ParsedExpr>, Vec<Witness>),
    OddSystem(OddSystemOp, Vec<ParsedExpr>, Vec<Witness>),
    MapTransfer(MapTransferOp, Vec<ParsedExpr>),
    Binary(BinaryOp, Box<ParsedExpr>, Box<ParsedExpr>),
    Lambda(Label, ParsedType, Box<ParsedExpr>),
    Call(Box<ParsedExpr>, Box<ParsedExpr>),
    Pair(Box<ParsedExpr>, Box<ParsedExpr>),
    Fst(Box<ParsedExpr>),
    Snd(Box<ParsedExpr>),
    Table(ParsedType, ParsedType, Entries),
    Restrict(Label, Box<ParsedExpr>),
    Sew(ParsedType, ParsedType, Entries),
}

pub(crate) struct Statement {
    pub(crate) span: Span,
    pub(crate) kind: StatementForm,
}

pub(crate) enum StatementForm {
    Space(Label, Vec<Vec<Label>>),
    QModule(Label, Vec<usize>, crate::q_module::Presentation),
    Let(Label, ParsedExpr),
    Show(ParsedExpr),
    Card(ParsedExpr),
    Assert(ParsedExpr, ParsedExpr),
    Betti(ParsedType),
    ExtDimension(usize, usize),
    PExtDimension(usize, usize, usize),
    PExtAnDimension(usize, usize, usize, usize),
    PMayDimension(usize, usize, usize, usize, usize),
    PMayAnDimension(usize, usize, usize, usize, usize, usize),
    ExtAnDimension(usize, usize, usize),
    XExtDimension(ParsedType, usize, usize),
    MapExtDimension(ParsedType, ParsedType, usize, i64),
    MapMayDimension(MapMaySpec),
    QExtDimension(QSpec),
    VExtDimension(QSpec<i64>),
    VBockDimension(VBockSpec),
    MayDimension(usize, usize, usize, usize),
    MayAnDimension(usize, usize, usize, usize, usize),
    CeDimension(usize, usize, usize, usize, usize),
    Homotopy(ParsedExpr, ParsedExpr),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Token<'src> {
    Word(&'src str),
    Number(&'src str),
    LBrace,
    RBrace,
    LParen,
    RParen,
    LBracket,
    RBracket,
    Colon,
    Semi,
    Comma,
    Dot,
    Less,
    Star,
    Slash,
    Equal,
    EqualEqual,
    Minus,
    Arrow,
    FatArrow,
    Tilde,
    Eof,
}

impl Token<'_> {
    fn description(self) -> String {
        match self {
            Self::Word(word) => format!("'{word}'"),
            Self::Number(number) => format!("'{number}'"),
            Self::Eof => "end of input".into(),
            symbol => {
                let spelling = match symbol {
                    Self::LBrace => "{",
                    Self::RBrace => "}",
                    Self::LParen => "(",
                    Self::RParen => ")",
                    Self::LBracket => "[",
                    Self::RBracket => "]",
                    Self::Colon => ":",
                    Self::Semi => ";",
                    Self::Comma => ",",
                    Self::Dot => ".",
                    Self::Less => "<",
                    Self::Star => "*",
                    Self::Slash => "\\",
                    Self::Equal => "=",
                    Self::EqualEqual => "==",
                    Self::Minus => "-",
                    Self::Arrow => "->",
                    Self::FatArrow => "=>",
                    Self::Tilde => "~",
                    Self::Word(_) | Self::Number(_) | Self::Eof => unreachable!(),
                };
                format!("'{spelling}'")
            }
        }
    }
}

struct Lexer<'src> {
    source: &'src str,
    offset: usize,
    span: Span,
}

impl<'src> Lexer<'src> {
    fn new(source: &'src str) -> Self {
        Self {
            source,
            offset: 0,
            span: Span { line: 1, column: 1 },
        }
    }

    fn peek(&self) -> Option<u8> {
        self.source.as_bytes().get(self.offset).copied()
    }

    fn step(&mut self) {
        if self.peek() == Some(b'\n') {
            self.span.line += 1;
            self.span.column = 1;
        } else {
            self.span.column += 1;
        }
        self.offset += 1;
    }

    fn eat(&mut self, byte: u8) -> bool {
        if self.peek() == Some(byte) {
            self.step();
            return true;
        }
        false
    }

    fn next(&mut self) -> Result<(Token<'src>, Span)> {
        loop {
            match self.peek() {
                Some(byte) if byte.is_ascii_whitespace() => self.step(),
                Some(b'/') if self.source.as_bytes().get(self.offset + 1) == Some(&b'/') => {
                    while self.peek().is_some_and(|byte| byte != b'\n') {
                        self.step();
                    }
                }
                _ => break,
            }
        }

        let span = self.span;
        let Some(byte) = self.peek() else {
            return Ok((Token::Eof, span));
        };
        let start = self.offset;
        self.step();
        let token = match byte {
            b'a'..=b'z' | b'A'..=b'Z' | b'_' => {
                while self
                    .peek()
                    .is_some_and(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
                {
                    self.step();
                }
                if self.offset - start > MAX_IDENTIFIER_BYTES {
                    return Err(Error::new(span, "identifier exceeds 128 bytes"));
                }
                Token::Word(&self.source[start..self.offset])
            }
            b'0'..=b'9' => {
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.step();
                }
                Token::Number(&self.source[start..self.offset])
            }
            b'{' => Token::LBrace,
            b'}' => Token::RBrace,
            b'(' => Token::LParen,
            b')' => Token::RParen,
            b'[' => Token::LBracket,
            b']' => Token::RBracket,
            b':' => Token::Colon,
            b';' => Token::Semi,
            b',' => Token::Comma,
            b'.' => Token::Dot,
            b'<' => Token::Less,
            b'*' => Token::Star,
            b'\\' => Token::Slash,
            b'=' if self.eat(b'=') => Token::EqualEqual,
            b'=' if self.eat(b'>') => Token::FatArrow,
            b'=' => Token::Equal,
            b'-' if self.eat(b'>') => Token::Arrow,
            b'-' => Token::Minus,
            b'~' => Token::Tilde,
            _ => {
                let character = self.source[start..].chars().next().unwrap();
                return Err(Error::new(
                    span,
                    format!("unexpected character {character:?}"),
                ));
            }
        };
        Ok((token, span))
    }
}

pub(crate) struct Parser<'src> {
    lexer: Lexer<'src>,
    lookahead: (Token<'src>, Span),
}

impl<'src> Parser<'src> {
    pub(crate) fn new(source: &'src str) -> Result<Self> {
        let mut lexer = Lexer::new(source);
        let lookahead = lexer.next()?;
        Ok(Self { lexer, lookahead })
    }

    fn token(&self) -> Token<'src> {
        self.lookahead.0
    }

    fn span(&self) -> Span {
        self.lookahead.1
    }

    fn advance(&mut self) -> Result<()> {
        self.lookahead = self.lexer.next()?;
        Ok(())
    }

    fn eat(&mut self, token: Token<'src>) -> Result<bool> {
        if self.token() == token {
            self.advance()?;
            return Ok(true);
        }
        Ok(false)
    }

    fn expect(&mut self, token: Token<'src>) -> Result<()> {
        if !self.eat(token)? {
            return Err(self.expected(&token.description()));
        }
        Ok(())
    }

    fn expected(&self, expected: &str) -> Error {
        Error::new(
            self.span(),
            format!("expected {expected}, found {}", self.token().description()),
        )
    }

    fn label(&mut self) -> Result<Label> {
        let Token::Word(text) = self.token() else {
            return Err(self.expected("an identifier"));
        };
        if matches!(
            text,
            "space"
                | "let"
                | "show"
                | "assert"
                | "betti"
                | "homotopy"
                | "Star"
                | "fst"
                | "snd"
                | "tab"
                | "res"
                | "sew"
        ) {
            return Err(self.expected("an identifier"));
        }
        let label = Label {
            text: text.into(),
            span: self.span(),
        };
        self.advance()?;
        Ok(label)
    }

    fn number(&mut self) -> Result<(usize, Span)> {
        let span = self.span();
        let Token::Number(text) = self.token() else {
            return Err(self.expected("a nonnegative decimal integer"));
        };
        let number = text
            .parse()
            .map_err(|_| Error::new(span, "decimal integer is out of range"))?;
        self.advance()?;
        Ok((number, span))
    }

    fn signed_degree(&mut self) -> Result<i64> {
        let span = self.span();
        let negative = self.eat(Token::Minus)?;
        let Token::Number(text) = self.token() else {
            return Err(self.expected("a signed internal degree"));
        };
        let out_of_range = || Error::new(span, "signed internal degree is out of range");
        let magnitude = i128::from(text.parse::<u64>().map_err(|_| out_of_range())?);
        let degree = i64::try_from(if negative { -magnitude } else { magnitude })
            .map_err(|_| out_of_range())?;
        self.advance()?;
        Ok(degree)
    }

    fn exponents(&mut self) -> Result<Vec<usize>> {
        self.expect(Token::LBracket)?;
        let mut exponents = Vec::new();
        while self.token() != Token::RBracket {
            exponents.push(self.number()?.0);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(exponents)
    }

    fn q_module_body(&mut self) -> Result<crate::q_module::Presentation> {
        self.expect(Token::LBrace)?;
        let mut generators = Vec::new();
        while self.token() != Token::RBrace {
            let name = self.label()?;
            self.expect(Token::Colon)?;
            generators.push((name, self.number()?.0));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::LBrace)?;
        let mut rules = Vec::new();
        while self.token() != Token::RBrace {
            let (primitive, span) = self.number()?;
            self.expect(Token::LParen)?;
            let source = self.label()?;
            self.expect(Token::RParen)?;
            self.expect(Token::FatArrow)?;
            self.expect(Token::LBracket)?;
            let mut targets = Vec::new();
            while self.token() != Token::RBracket {
                targets.push(self.label()?);
                if !self.eat(Token::Comma)? {
                    break;
                }
            }
            self.expect(Token::RBracket)?;
            rules.push(crate::q_module::Rule {
                primitive,
                span,
                source,
                targets,
            });
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(crate::q_module::Presentation { generators, rules })
    }

    fn chains(&mut self) -> Result<Vec<Vec<Label>>> {
        self.expect(Token::LBrace)?;
        let mut chains = Vec::new();
        while self.token() != Token::RBrace {
            let mut chain = vec![self.label()?];
            while self.eat(Token::Less)? {
                chain.push(self.label()?);
            }
            chains.push(chain);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(chains)
    }

    pub(crate) fn statement(&mut self) -> Result<Option<Statement>> {
        let span = self.span();
        let kind = match self.token() {
            Token::Eof => return Ok(None),
            Token::Word("qmodule") => {
                self.advance()?;
                let name = self.label()?;
                let selected = self.exponents()?;
                let presentation = self.q_module_body()?;
                StatementForm::QModule(name, selected, presentation)
            }
            Token::Word("space") => {
                self.advance()?;
                let name = self.label()?;
                StatementForm::Space(name, self.chains()?)
            }
            Token::Word("let") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::Equal)?;
                StatementForm::Let(name, self.expr(0)?)
            }
            Token::Word("show") => {
                self.advance()?;
                StatementForm::Show(self.expr(0)?)
            }
            Token::Word("card") => {
                self.advance()?;
                StatementForm::Card(self.expr(0)?)
            }
            Token::Word("assert") => {
                self.advance()?;
                let left = self.expr(0)?;
                self.expect(Token::EqualEqual)?;
                StatementForm::Assert(left, self.expr(0)?)
            }
            Token::Word("betti") => {
                self.advance()?;
                StatementForm::Betti(self.ty(0)?)
            }
            Token::Word("pextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (p, s, t) = self.prime_bidegree()?;
                self.expect(Token::RBracket)?;
                StatementForm::PExtDimension(p, s, t)
            }
            Token::Word("pextandim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (p, n, s, t) = self.prime_profile_bidegree()?;
                self.expect(Token::RBracket)?;
                StatementForm::PExtAnDimension(p, n, s, t)
            }
            Token::Word("extdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (s, t) = self.bidegree()?;
                self.expect(Token::RBracket)?;
                StatementForm::ExtDimension(s, t)
            }
            Token::Word("extandim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (n, s, t) = self.profile_bidegree()?;
                self.expect(Token::RBracket)?;
                StatementForm::ExtAnDimension(n, s, t)
            }
            Token::Word("xextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (base, s, t) = self.space_bidegree(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::XExtDimension(base, s, t)
            }
            Token::Word("mapextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (source, target, s, t) = self.map_bidegree(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::MapExtDimension(source, target, s, t)
            }
            Token::Word("qextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.q_spec(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::QExtDimension(spec)
            }
            Token::Word("vextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.v_spec(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::VExtDimension(spec)
            }
            Token::Word("vbockdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.vbock_spec(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::VBockDimension(spec)
            }
            Token::Word("mapmaydim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.map_may_spec(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::MapMayDimension(spec)
            }
            Token::Word(name @ ("mayandim" | "ceandim")) => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (n, r, s, t, f) = self.profile_page_degree()?;
                self.expect(Token::RBracket)?;
                if name == "ceandim" {
                    StatementForm::CeDimension(n, r, s, t, f)
                } else {
                    StatementForm::MayAnDimension(n, r, s, t, f)
                }
            }
            Token::Word("maydim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (r, s, t, f) = self.may_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::MayDimension(r, s, t, f)
            }
            Token::Word("pmaydim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (p, r, s, t, f) = self.pmay_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::PMayDimension(p, r, s, t, f)
            }
            Token::Word("pmayandim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (p, n, r, s, t, f) = self.prime_profile_page_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::PMayAnDimension(p, n, r, s, t, f)
            }
            Token::Word("homotopy") => {
                self.advance()?;
                let left = self.expr(0)?;
                self.expect(Token::Tilde)?;
                StatementForm::Homotopy(left, self.expr(0)?)
            }
            _ => return Err(self.expected("a statement")),
        };
        self.expect(Token::Semi)?;
        Ok(Some(Statement { span, kind }))
    }

    fn prime_bidegree(&mut self) -> Result<(usize, usize, usize)> {
        let (p, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (s, t) = self.bidegree()?;
        Ok((p, s, t))
    }

    fn prime_profile_bidegree(&mut self) -> Result<(usize, usize, usize, usize)> {
        let (p, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (n, s, t) = self.profile_bidegree()?;
        Ok((p, n, s, t))
    }

    fn bidegree(&mut self) -> Result<(usize, usize)> {
        let (s, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (t, _) = self.number()?;
        Ok((s, t))
    }

    fn profile_bidegree(&mut self) -> Result<(usize, usize, usize)> {
        let (n, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (s, t) = self.bidegree()?;
        Ok((n, s, t))
    }

    fn space_bidegree(&mut self, nest: usize) -> Result<(ParsedType, usize, usize)> {
        let base = self.ty(nest + 1)?;
        self.expect(Token::Comma)?;
        let (s, t) = self.bidegree()?;
        Ok((base, s, t))
    }

    fn map_endpoints(&mut self, nest: usize) -> Result<(ParsedType, ParsedType)> {
        let source = self.ty(nest + 1)?;
        self.expect(Token::Comma)?;
        let target = self.ty(nest + 1)?;
        self.expect(Token::Comma)?;
        Ok((source, target))
    }

    fn map_bidegree(&mut self, nest: usize) -> Result<(ParsedType, ParsedType, usize, i64)> {
        let (source, target) = self.map_endpoints(nest)?;
        let (s, _) = self.number()?;
        self.expect(Token::Comma)?;
        let t = self.signed_degree()?;
        Ok((source, target, s, t))
    }

    fn q_spec(&mut self, nest: usize) -> Result<QSpec> {
        self.koszul_spec(nest, |parser| parser.number().map(|(s, _)| s))
    }

    fn v_spec(&mut self, nest: usize) -> Result<QSpec<i64>> {
        self.koszul_spec(nest, Self::signed_degree)
    }

    fn vbock_spec(&mut self, nest: usize) -> Result<VBockSpec> {
        let QSpec {
            source,
            target,
            selected,
            s: (r, s),
            t,
        } = self.koszul_spec(nest, |parser| {
            let (r, _) = parser.number()?;
            parser.expect(Token::Comma)?;
            Ok((r, parser.signed_degree()?))
        })?;
        self.expect(Token::Comma)?;
        let (f, _) = self.number()?;
        Ok(VBockSpec {
            localized: QSpec {
                source,
                target,
                selected,
                s,
                t,
            },
            r,
            f,
        })
    }

    fn koszul_spec<S>(
        &mut self,
        nest: usize,
        filtration: impl FnOnce(&mut Self) -> Result<S>,
    ) -> Result<QSpec<S>> {
        let (source, target) = self.map_endpoints(nest)?;
        let selected = self.exponents()?;
        self.expect(Token::Comma)?;
        let s = filtration(self)?;
        self.expect(Token::Comma)?;
        let t = self.signed_degree()?;
        Ok(QSpec {
            source: Box::new(source),
            target: Box::new(target),
            selected,
            s,
            t,
        })
    }

    fn map_may_spec(&mut self, nest: usize) -> Result<MapMaySpec> {
        let (source, target) = self.map_endpoints(nest)?;
        let (r, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (s, _) = self.number()?;
        self.expect(Token::Comma)?;
        let t = self.signed_degree()?;
        self.expect(Token::Comma)?;
        let (f, _) = self.number()?;
        Ok(MapMaySpec {
            source: Box::new(source),
            target: Box::new(target),
            r,
            s,
            t,
            f,
        })
    }

    fn map_generator(&mut self, nest: usize) -> Result<MapGenerator> {
        let base = self.ty(nest + 1)?;
        self.expect(Token::Comma)?;
        let (degree, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (index, span) = self.number()?;
        Ok(MapGenerator {
            base,
            degree,
            index,
            span,
        })
    }

    fn profile_page_degree(&mut self) -> Result<(usize, usize, usize, usize, usize)> {
        let (n, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (r, s, t, f) = self.may_degree()?;
        Ok((n, r, s, t, f))
    }

    fn may_degree(&mut self) -> Result<(usize, usize, usize, usize)> {
        let (r, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (s, t) = self.bidegree()?;
        self.expect(Token::Comma)?;
        let (f, _) = self.number()?;
        Ok((r, s, t, f))
    }

    fn pmay_degree(&mut self) -> Result<(usize, usize, usize, usize, usize)> {
        let (prime, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (r, s, t, f) = self.may_degree()?;
        Ok((prime, r, s, t, f))
    }

    fn prime_profile_page_degree(&mut self) -> Result<(usize, usize, usize, usize, usize, usize)> {
        let (prime, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (n, r, s, t, f) = self.profile_page_degree()?;
        Ok((prime, n, r, s, t, f))
    }

    fn ty(&mut self, nest: usize) -> Result<ParsedType> {
        check_nesting(nest, self.span())?;
        let mut head = self.type_atom(nest)?;
        while self.eat(Token::Star)? {
            let right = self.type_atom(nest)?;
            let span = head.span;
            let depth = head.depth.max(right.depth) + 1;
            head = ty(
                span,
                TypeForm::Product(Box::new(head), Box::new(right)),
                depth,
            )?;
        }
        if self.eat(Token::Arrow)? {
            let right = self.ty(nest + 1)?;
            let span = head.span;
            let depth = head.depth.max(right.depth) + 1;
            head = ty(
                span,
                TypeForm::Arrow(Box::new(head), Box::new(right)),
                depth,
            )?;
        }
        Ok(head)
    }

    fn type_atom(&mut self, nest: usize) -> Result<ParsedType> {
        let span = self.span();
        if self.eat(Token::LParen)? {
            let ty = self.ty(nest + 1)?;
            self.expect(Token::RParen)?;
            return Ok(ty);
        }
        if self.eat(Token::Word("Star"))? {
            self.expect(Token::LParen)?;
            let base = self.ty(nest + 1)?;
            self.expect(Token::Comma)?;
            let point = self.label()?;
            self.expect(Token::RParen)?;
            let depth = base.depth + 1;
            return ty(span, TypeForm::Star(Box::new(base), point), depth);
        }
        let name = self.label()?;
        if matches!(name.text.as_str(), "PCobarAn" | "PExtAn") && self.eat(Token::LParen)? {
            let (p, n, s, t) = self.prime_profile_bidegree()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "PCobarAn" {
                TypeForm::PCobarAn(p, n, s, t)
            } else {
                TypeForm::PExtAn(p, n, s, t)
            };
            return ty(span, kind, 0);
        }
        if matches!(name.text.as_str(), "PCobar" | "PExt") && self.eat(Token::LParen)? {
            let (p, s, t) = self.prime_bidegree()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "PCobar" {
                TypeForm::PCobar(p, s, t)
            } else {
                TypeForm::PExt(p, s, t)
            };
            return ty(span, kind, 0);
        }
        if matches!(name.text.as_str(), "CobarAn" | "ExtAn") && self.eat(Token::LParen)? {
            let (n, s, t) = self.profile_bidegree()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "CobarAn" {
                TypeForm::CobarAn(n, s, t)
            } else {
                TypeForm::ExtAn(n, s, t)
            };
            return ty(span, kind, 0);
        }
        if name.text == "XExt" && self.eat(Token::LParen)? {
            let (base, s, t) = self.space_bidegree(nest)?;
            self.expect(Token::RParen)?;
            let depth = base.depth + 1;
            return ty(span, TypeForm::XExt(Box::new(base), s, t), depth);
        }
        if matches!(name.text.as_str(), "MapCobar" | "MapExt") && self.eat(Token::LParen)? {
            let (source, target, s, t) = self.map_bidegree(nest)?;
            self.expect(Token::RParen)?;
            let depth = source.depth.max(target.depth) + 1;
            let source = Box::new(source);
            let target = Box::new(target);
            let kind = if name.text == "MapCobar" {
                TypeForm::MapCobar(source, target, s, t)
            } else {
                TypeForm::MapExt(source, target, s, t)
            };
            return ty(span, kind, depth);
        }
        if matches!(name.text.as_str(), "QKoszul" | "QExt") && self.eat(Token::LParen)? {
            let spec = self.q_spec(nest)?;
            self.expect(Token::RParen)?;
            let depth = spec.source.depth.max(spec.target.depth) + 1;
            let kind = if name.text == "QKoszul" {
                TypeForm::QKoszul(spec)
            } else {
                TypeForm::QExt(spec)
            };
            return ty(span, kind, depth);
        }
        if matches!(name.text.as_str(), "VKoszul" | "VExt") && self.eat(Token::LParen)? {
            let spec = self.v_spec(nest)?;
            self.expect(Token::RParen)?;
            let depth = spec.source.depth.max(spec.target.depth) + 1;
            let kind = if name.text == "VKoszul" {
                TypeForm::VKoszul(spec)
            } else {
                TypeForm::VExt(spec)
            };
            return ty(span, kind, depth);
        }
        if name.text == "VBock" && self.eat(Token::LParen)? {
            let spec = self.vbock_spec(nest)?;
            self.expect(Token::RParen)?;
            let depth = spec.localized.source.depth.max(spec.localized.target.depth) + 1;
            return ty(span, TypeForm::VBock(spec), depth);
        }
        if name.text == "MapMay" && self.eat(Token::LParen)? {
            let spec = self.map_may_spec(nest)?;
            self.expect(Token::RParen)?;
            let depth = spec.source.depth.max(spec.target.depth) + 1;
            return ty(span, TypeForm::MapMay(spec), depth);
        }
        if matches!(name.text.as_str(), "MayAn" | "CEAn") && self.eat(Token::LParen)? {
            let (n, r, s, t, f) = self.profile_page_degree()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "CEAn" {
                TypeForm::CeAn(n, r, s, t, f)
            } else {
                TypeForm::MayAn(n, r, s, t, f)
            };
            return ty(span, kind, 0);
        }
        if name.text == "May" && self.eat(Token::LParen)? {
            let (r, s, t, f) = self.may_degree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::May(r, s, t, f), 0);
        }
        if name.text == "PMay" && self.eat(Token::LParen)? {
            let (p, r, s, t, f) = self.pmay_degree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::PMay(p, r, s, t, f), 0);
        }
        if name.text == "PMayAn" && self.eat(Token::LParen)? {
            let (p, n, r, s, t, f) = self.prime_profile_page_degree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::PMayAn(p, n, r, s, t, f), 0);
        }
        if matches!(
            name.text.as_str(),
            "A" | "Bar" | "H2" | "H" | "Susp" | "Stable" | "Cobar" | "Ext"
        ) && self.eat(Token::LParen)?
        {
            let (kind, depth) = if matches!(name.text.as_str(), "A" | "Bar" | "H2") {
                let (degree, _) = self.number()?;
                let kind = match name.text.as_str() {
                    "Bar" => TypeForm::Bar(degree),
                    "H2" => TypeForm::CyclicH2(degree),
                    _ => TypeForm::Operations(degree),
                };
                (kind, 0)
            } else if matches!(name.text.as_str(), "Cobar" | "Ext") {
                let (s, t) = self.bidegree()?;
                let kind = match name.text.as_str() {
                    "Cobar" => TypeForm::Cobar(s, t),
                    _ => TypeForm::Ext(s, t),
                };
                (kind, 0)
            } else {
                let base = self.ty(nest + 1)?;
                let depth = base.depth + 1;
                let kind = if name.text == "Susp" {
                    TypeForm::Suspension(Box::new(base))
                } else {
                    self.expect(Token::Comma)?;
                    let (first, _) = self.number()?;
                    if name.text == "H" {
                        TypeForm::Cohomology(Box::new(base), first)
                    } else {
                        self.expect(Token::Comma)?;
                        let (degree, _) = self.number()?;
                        TypeForm::Stable(Box::new(base), first, degree)
                    }
                };
                (kind, depth)
            };
            self.expect(Token::RParen)?;
            return ty(span, kind, depth);
        }
        ty(span, TypeForm::Name(name), 0)
    }

    fn expr(&mut self, nest: usize) -> Result<ParsedExpr> {
        check_nesting(nest, self.span())?;
        let mut head = self.expr_atom(nest)?;
        while self.eat(Token::LParen)? {
            if let ExprForm::Name(name) = &head.kind {
                if matches!(name.text.as_str(), "massey" | "powerdef" | "bracket")
                    || TransferOp::from_name(&name.text).is_some()
                    || MapBracketOp::from_name(&name.text).is_some()
                    || MapSystemOp::from_name(&name.text).is_some()
                    || OddSystemOp::from_name(&name.text).is_some()
                    || MapTransferOp::from_name(&name.text).is_some()
                {
                    head = self.variadic(head, nest)?;
                    continue;
                }
            }
            let arg = self.expr(nest + 1)?;
            if self.eat(Token::Comma)? {
                let primitive = match &head.kind {
                    ExprForm::Name(name) => BinaryOp::from_name(&name.text),
                    _ => None,
                };
                let Some(primitive) = primitive else {
                    return Err(Error::new(
                        head.span,
                        "a function call takes one argument; use nested calls or pass a pair",
                    ));
                };
                let right = self.expr(nest + 1)?;
                self.expect(Token::RParen)?;
                let depth = arg.depth.max(right.depth) + 1;
                head = expr(
                    head.span,
                    ExprForm::Binary(primitive, Box::new(arg), Box::new(right)),
                    depth,
                )?;
                continue;
            }
            self.expect(Token::RParen)?;
            let span = head.span;
            let depth = head.depth.max(arg.depth) + 1;
            head = expr(span, ExprForm::Call(Box::new(head), Box::new(arg)), depth)?;
        }
        Ok(head)
    }

    fn arguments(&mut self, nest: usize) -> Result<Vec<ParsedExpr>> {
        let mut arguments = Vec::new();
        while self.token() != Token::RParen {
            arguments.push(self.expr(nest + 1)?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RParen)?;
        Ok(arguments)
    }

    fn variadic(&mut self, head: ParsedExpr, nest: usize) -> Result<ParsedExpr> {
        let mut arguments = self.arguments(nest)?;
        let ExprForm::Name(name) = &head.kind else {
            unreachable!();
        };
        let massey = name.text == "massey";
        let map_system = MapSystemOp::from_name(&name.text);
        let odd_system = OddSystemOp::from_name(&name.text);
        let has_witnesses =
            massey || map_system.is_some() || odd_system.is_some_and(OddSystemOp::explicit);
        let span = head.span;
        let mut depth = arguments
            .iter()
            .map(|argument| argument.depth + 1)
            .max()
            .unwrap_or(0);
        if arguments.len() == 1 && (self.token() != Token::LBrace || !has_witnesses) {
            return expr(
                span,
                ExprForm::Call(Box::new(head), Box::new(arguments.pop().unwrap())),
                depth,
            );
        }
        if let Some(operation) = TransferOp::from_name(&name.text) {
            return expr(span, ExprForm::Transfer(operation, arguments), depth);
        }
        if let Some(operation) = MapBracketOp::from_name(&name.text) {
            return expr(span, ExprForm::MapBracket(operation, arguments), depth);
        }
        if let Some(operation) = MapTransferOp::from_name(&name.text) {
            return expr(span, ExprForm::MapTransfer(operation, arguments), depth);
        }
        if let Some(operation) = map_system {
            let witnesses = self.witnesses(nest, &mut depth)?;
            return expr(
                span,
                ExprForm::MapSystem(operation, arguments, witnesses),
                depth,
            );
        }
        if let Some(operation) = odd_system {
            let witnesses = if operation.explicit() {
                self.witnesses(nest, &mut depth)?
            } else {
                Vec::new()
            };
            return expr(
                span,
                ExprForm::OddSystem(operation, arguments, witnesses),
                depth,
            );
        }
        if name.text == "bracket" {
            return expr(span, ExprForm::Bracket(arguments), depth);
        }
        if !massey {
            return expr(span, ExprForm::PowerDef(arguments), depth);
        }

        let witnesses = self.witnesses(nest, &mut depth)?;
        expr(span, ExprForm::Massey(arguments, witnesses), depth)
    }

    fn witnesses(&mut self, nest: usize, depth: &mut usize) -> Result<Vec<Witness>> {
        self.expect(Token::LBrace)?;
        let mut witnesses = Vec::new();
        while self.token() != Token::RBrace {
            let (from, span) = self.number()?;
            self.expect(Token::Less)?;
            let (to, _) = self.number()?;
            self.expect(Token::FatArrow)?;
            let value = self.expr(nest + 1)?;
            *depth = (*depth).max(value.depth + 1);
            witnesses.push(Witness {
                from,
                to,
                span,
                value,
            });
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(witnesses)
    }

    fn expr_atom(&mut self, nest: usize) -> Result<ParsedExpr> {
        let span = self.span();
        match self.token() {
            Token::Slash => {
                self.advance()?;
                let parameter = self.label()?;
                self.expect(Token::Colon)?;
                let domain = self.ty(nest + 1)?;
                self.expect(Token::FatArrow)?;
                let body = self.expr(nest + 1)?;
                let depth = body.depth + 1;
                expr(
                    span,
                    ExprForm::Lambda(parameter, domain, Box::new(body)),
                    depth,
                )
            }
            Token::LParen => {
                self.advance()?;
                let left = self.expr(nest + 1)?;
                if !self.eat(Token::Comma)? {
                    self.expect(Token::RParen)?;
                    return Ok(left);
                }
                let right = self.expr(nest + 1)?;
                self.expect(Token::RParen)?;
                let depth = left.depth.max(right.depth) + 1;
                expr(span, ExprForm::Pair(Box::new(left), Box::new(right)), depth)
            }
            Token::Word("fst") | Token::Word("snd") => {
                let first = self.token() == Token::Word("fst");
                self.advance()?;
                let value = self.paren_expr(nest)?;
                let depth = value.depth + 1;
                let kind = if first {
                    ExprForm::Fst(Box::new(value))
                } else {
                    ExprForm::Snd(Box::new(value))
                };
                expr(span, kind, depth)
            }
            Token::Word("res") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let point = self.label()?;
                self.expect(Token::RBracket)?;
                let base = self.paren_expr(nest)?;
                let depth = base.depth + 1;
                expr(span, ExprForm::Restrict(point, Box::new(base)), depth)
            }
            Token::Word("tab") | Token::Word("sew") => {
                let is_table = self.token() == Token::Word("tab");
                let keyword = self.token().description();
                self.advance()?;
                let ty = self.ty(nest + 1)?;
                let TypeForm::Arrow(domain, codomain) = ty.kind else {
                    return Err(Error::new(
                        ty.span,
                        format!("{keyword} requires a type of the form X -> Y"),
                    ));
                };
                self.expect(Token::LBrace)?;
                let mut entries = Vec::new();
                while self.token() != Token::RBrace {
                    let label = self.label()?;
                    self.expect(Token::FatArrow)?;
                    entries.push((label, self.expr(nest + 1)?));
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                let depth = entries
                    .iter()
                    .map(|(_, value)| value.depth + 1)
                    .max()
                    .unwrap_or(0);
                let kind = if is_table {
                    ExprForm::Table(*domain, *codomain, entries)
                } else {
                    ExprForm::Sew(*domain, *codomain, entries)
                };
                expr(span, kind, depth)
            }
            Token::Word("Star") => {
                let ty = self.type_atom(nest)?;
                self.expect(Token::Dot)?;
                expr(span, ExprForm::Point(ty, self.label()?), 0)
            }
            Token::Word(_) => {
                let name = self.label()?;
                if name.text == "milnor" && self.eat(Token::LBrace)? {
                    let exponents = self.exponents()?;
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::Milnor(exponents), 0);
                }
                if name.text == "cword" && self.eat(Token::LBrace)? {
                    let mut factors = Vec::new();
                    while self.token() != Token::RBrace {
                        let at = self.span();
                        factors.push((self.exponents()?, at));
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::CobarWord(factors), 0);
                }
                if name.text == "pword" && self.eat(Token::LBracket)? {
                    let (p, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut factors = Vec::new();
                    while self.token() != Token::RBrace {
                        let at = self.span();
                        self.expect(Token::LParen)?;
                        let exponents = self.exponents()?;
                        self.expect(Token::Comma)?;
                        let taus = self.exponents()?;
                        self.expect(Token::RParen)?;
                        factors.push((exponents, taus, at));
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::PCobarWord(p, factors), 0);
                }
                if name.text == "pscale" && self.eat(Token::LBracket)? {
                    let (scalar, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::PScale(scalar, Box::new(value)), depth);
                }
                if name.text == "pextgen" && self.eat(Token::LBracket)? {
                    let (p, s, t) = self.prime_bidegree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::PExtGenerator(p, s, t, index, at), 0);
                }
                if name.text == "pextangen" && self.eat(Token::LBracket)? {
                    let (p, n, s, t) = self.prime_profile_bidegree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::PExtAnGenerator(p, n, s, t, index, at), 0);
                }
                if name.text == "extgen" && self.eat(Token::LBracket)? {
                    let (s, t) = self.bidegree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::ExtGenerator(s, t, index, at), 0);
                }
                if name.text == "extangen" && self.eat(Token::LBracket)? {
                    let (n, s, t) = self.profile_bidegree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::ExtAnGenerator(n, s, t, index, at), 0);
                }
                if name.text == "resA" && self.eat(Token::LBracket)? {
                    let (n, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::RestrictAlgebra(n, Box::new(value)), depth);
                }
                if name.text == "cCup" && self.eat(Token::LBracket)? {
                    let (i, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LParen)?;
                    let left = self.expr(nest + 1)?;
                    self.expect(Token::Comma)?;
                    let right = self.expr(nest + 1)?;
                    self.expect(Token::RParen)?;
                    let depth = left.depth.max(right.depth) + 1;
                    return expr(
                        span,
                        ExprForm::CobarCup(i, Box::new(left), Box::new(right)),
                        depth,
                    );
                }
                if name.text == "extSq" && self.eat(Token::LBracket)? {
                    let (k, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::ExtSquare(k, Box::new(value)), depth);
                }
                if name.text == "cNakamura" && self.eat(Token::LBracket)? {
                    let (k, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::Nakamura(k, Box::new(value)), depth);
                }
                if name.text == "mayLiftTo" && self.eat(Token::LBracket)? {
                    let (r, f) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::MayLiftTo(r, f, Box::new(value)), depth);
                }
                if name.text == "resHigher" && self.eat(Token::LBracket)? {
                    let (n, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LParen)?;
                    let inputs = self.arguments(nest)?;
                    let depth = inputs
                        .iter()
                        .map(|input| input.depth + 1)
                        .max()
                        .unwrap_or(0);
                    return expr(span, ExprForm::HigherRestriction(n, inputs), depth);
                }
                if name.text == "xextgen" && self.eat(Token::LBracket)? {
                    let (base, s, t) = self.space_bidegree(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::XExtGenerator(base, s, t, index, at), 0);
                }
                if name.text == "mapextgen" && self.eat(Token::LBracket)? {
                    let (source, target, s, t) = self.map_bidegree(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::MapExtGenerator(source, target, s, t, index, at),
                        0,
                    );
                }
                if name.text == "mword" && self.eat(Token::LBracket)? {
                    let source = self.map_generator(nest)?;
                    self.expect(Token::Comma)?;
                    let target = self.map_generator(nest)?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::MapWord(source, target, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "eid" && self.eat(Token::LBracket)? {
                    let base = self.ty(nest + 1)?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MapExtIdentity(base), 0);
                }
                if name.text == "qterm" && self.eat(Token::LBracket)? {
                    let source = self.map_generator(nest)?;
                    self.expect(Token::Comma)?;
                    let target = self.map_generator(nest)?;
                    self.expect(Token::Comma)?;
                    let selected = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut powers = Vec::new();
                    while self.token() != Token::RBrace {
                        powers.push(self.number()?.0);
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(
                        span,
                        ExprForm::QTerm(QTerm {
                            source,
                            target,
                            selected,
                            powers,
                        }),
                        0,
                    );
                }
                if name.text == "qid" && self.eat(Token::LBracket)? {
                    let base = self.ty(nest + 1)?;
                    self.expect(Token::Comma)?;
                    let selected = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::QIdentity(base, selected), 0);
                }
                if name.text == "qextgen" && self.eat(Token::LBracket)? {
                    let spec = self.q_spec(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::QExtGenerator(spec, index, at), 0);
                }
                if name.text == "vextgen" && self.eat(Token::LBracket)? {
                    let spec = self.v_spec(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::VExtGenerator(spec, index, at), 0);
                }
                if name.text == "vbockgen" && self.eat(Token::LBracket)? {
                    let spec = self.vbock_spec(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::VBockGenerator(spec, index, at), 0);
                }
                if name.text == "vbClass" && self.eat(Token::LBracket)? {
                    let (r, f) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::VBockClass(r, f, Box::new(value)), depth);
                }
                if name.text == "vShift" && self.eat(Token::LBracket)? {
                    let k = self.signed_degree()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::VShift(k, Box::new(input)), depth);
                }
                if name.text == "qRes" && self.token() == Token::LBracket {
                    let selected = self.exponents()?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::QRestrict(selected, Box::new(input)), depth);
                }
                if name.text == "mapmaygen" && self.eat(Token::LBracket)? {
                    let spec = self.map_may_spec(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MapMayGenerator(spec, index, at), 0);
                }
                if matches!(name.text.as_str(), "mayangen" | "ceangen")
                    && self.eat(Token::LBracket)?
                {
                    let (n, r, s, t, f) = self.profile_page_degree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let kind = if name.text == "ceangen" {
                        ExprForm::CeGenerator(n, r, s, t, f, index, at)
                    } else {
                        ExprForm::MayAnGenerator(n, r, s, t, f, index, at)
                    };
                    return expr(span, kind, 0);
                }
                if name.text == "maygen" && self.eat(Token::LBracket)? {
                    let (r, s, t, f) = self.may_degree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MayGenerator(r, s, t, f, index, at), 0);
                }
                if name.text == "pmaygen" && self.eat(Token::LBracket)? {
                    let (p, r, s, t, f) = self.pmay_degree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::PMayGenerator(p, r, s, t, f, index, at), 0);
                }
                if name.text == "pmayangen" && self.eat(Token::LBracket)? {
                    let (p, n, r, s, t, f) = self.prime_profile_page_degree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::PMayAnGenerator(p, n, r, s, t, f, index, at),
                        0,
                    );
                }
                if matches!(name.text.as_str(), "mayClass" | "ceClass")
                    && self.eat(Token::LBracket)?
                {
                    let (r, f) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    let kind = if name.text == "ceClass" {
                        ExprForm::CeClass(r, f, Box::new(value))
                    } else {
                        ExprForm::MayClass(r, f, Box::new(value))
                    };
                    return expr(span, kind, depth);
                }
                if matches!(name.text.as_str(), "bar" | "binom" | "carry")
                    && self.eat(Token::LBracket)?
                {
                    let (modulus, at) = self.number()?;
                    let form = match name.text.as_str() {
                        "binom" => {
                            self.expect(Token::Comma)?;
                            let (index, index_at) = self.number()?;
                            self.expect(Token::RBracket)?;
                            CyclicForm::Binomial(index, index_at)
                        }
                        "bar" => {
                            self.expect(Token::RBracket)?;
                            self.expect(Token::LBrace)?;
                            let mut support = Vec::new();
                            while self.token() != Token::RBrace {
                                support.push(self.number()?);
                                if !self.eat(Token::Comma)? {
                                    break;
                                }
                            }
                            self.expect(Token::RBrace)?;
                            CyclicForm::Support(support)
                        }
                        _ => {
                            self.expect(Token::RBracket)?;
                            CyclicForm::Carry
                        }
                    };
                    return expr(span, ExprForm::Cyclic(modulus, at, form), 0);
                }
                if name.text == "sq" && self.eat(Token::LBracket)? {
                    let (degree, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::Square(degree), 0);
                }
                if name.text == "coact" && self.eat(Token::LBracket)? {
                    let (n, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::Coact(n, index, at, Box::new(value)), depth);
                }
                if name.text == "shift" && self.eat(Token::LBracket)? {
                    let (shift, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::Shift(shift, Box::new(value)), depth);
                }
                let suspension = match name.text.as_str() {
                    "north" => Some(SuspensionForm::North),
                    "south" => Some(SuspensionForm::South),
                    "equator" => Some(SuspensionForm::Equator),
                    _ => None,
                };
                if let Some(form) = suspension {
                    if self.eat(Token::LBracket)? {
                        let base = self.ty(nest + 1)?;
                        self.expect(Token::RBracket)?;
                        return expr(span, ExprForm::Suspension(base, form), 0);
                    }
                }
                if matches!(name.text.as_str(), "zero" | "gen" | "cocycle")
                    && self.eat(Token::LBracket)?
                {
                    let base = self.ty(nest + 1)?;
                    if name.text == "zero" && self.eat(Token::RBracket)? {
                        return expr(span, ExprForm::Zero(base), 0);
                    }
                    self.expect(Token::Comma)?;
                    let (degree, _) = self.number()?;
                    let form = match name.text.as_str() {
                        "gen" => {
                            self.expect(Token::Comma)?;
                            let (index, at) = self.number()?;
                            self.expect(Token::RBracket)?;
                            ClassForm::Generator(index, at)
                        }
                        "cocycle" => {
                            self.expect(Token::RBracket)?;
                            ClassForm::Cocycle(self.chains()?)
                        }
                        _ => {
                            self.expect(Token::RBracket)?;
                            ClassForm::Zero
                        }
                    };
                    return expr(span, ExprForm::Class(base, degree, form), 0);
                }
                if self.eat(Token::Dot)? {
                    let ty = ty(span, TypeForm::Name(name), 0)?;
                    return expr(span, ExprForm::Point(ty, self.label()?), 0);
                }
                expr(span, ExprForm::Name(name), 0)
            }
            _ => Err(self.expected("an expression")),
        }
    }

    fn paren_expr(&mut self, nest: usize) -> Result<ParsedExpr> {
        self.expect(Token::LParen)?;
        let value = self.expr(nest + 1)?;
        self.expect(Token::RParen)?;
        Ok(value)
    }
}

fn check_nesting(depth: usize, span: Span) -> Result<()> {
    if depth > MAX_NESTING {
        return Err(Error::new(span, "syntax nesting exceeds 64 levels"));
    }
    Ok(())
}

fn ty(span: Span, kind: TypeForm, depth: usize) -> Result<ParsedType> {
    check_nesting(depth, span)?;
    Ok(ParsedType { span, kind, depth })
}

fn expr(span: Span, kind: ExprForm, depth: usize) -> Result<ParsedExpr> {
    check_nesting(depth, span)?;
    Ok(ParsedExpr { span, kind, depth })
}

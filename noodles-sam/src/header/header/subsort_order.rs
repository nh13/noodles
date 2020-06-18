use std::{error, fmt, str::FromStr};

/// A SAM header header subsort order (`SS`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SubsortOrder {
    /// Alignments are primarily unsorted (`unsorted`).
    Unsorted(String),
    /// Alignments are primarily sorted by read name (`queryname`).
    QueryName(String),
    /// Alignments are primarily sorted by reference sequence and position (`coordinate`).
    Coordinate(String),
}

impl fmt::Display for SubsortOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unsorted(subsort) => write!(f, "unsorted:{}", subsort),
            Self::QueryName(subsort) => write!(f, "queryname:{}", subsort),
            Self::Coordinate(subsort) => write!(f, "coordinate:{}", subsort),
        }
    }
}

/// An error returned when a raw SAM header header subsort order fails to parse.
#[derive(Debug)]
pub enum ParseError {
    /// The primary sort order is missing.
    MissingOrder,
    /// The primary sort order is invalid.
    InvalidOrder(String),
    /// The subsort order is missing.
    MissingSubsort,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid subsort order: ")?;

        match self {
            Self::MissingOrder => write!(f, "missing order"),
            Self::InvalidOrder(s) => {
                write!(f, "expected {{unsorted, queryname, coordinate}}, got {}", s)
            }
            Self::MissingSubsort => write!(f, "missing subsort"),
        }
    }
}

impl FromStr for SubsortOrder {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut pieces = s.splitn(2, ':');

        let order = pieces.next().ok_or_else(|| ParseError::MissingOrder)?;

        let subsort = pieces
            .next()
            .map(|s| s.into())
            .ok_or_else(|| ParseError::MissingSubsort)?;

        match order {
            "unsorted" => Ok(Self::Unsorted(subsort)),
            "queryname" => Ok(Self::QueryName(subsort)),
            "coordinate" => Ok(Self::Coordinate(subsort)),
            _ => Err(ParseError::InvalidOrder(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(
            SubsortOrder::Unsorted(String::from("MI")).to_string(),
            "unsorted:MI"
        );

        assert_eq!(
            SubsortOrder::QueryName(String::from("MI")).to_string(),
            "queryname:MI"
        );

        assert_eq!(
            SubsortOrder::Coordinate(String::from("MI")).to_string(),
            "coordinate:MI"
        );
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!(
            "unsorted:MI".parse::<SubsortOrder>()?,
            SubsortOrder::Unsorted(String::from("MI"))
        );

        assert_eq!(
            "queryname:MI".parse::<SubsortOrder>()?,
            SubsortOrder::QueryName(String::from("MI"))
        );

        assert_eq!(
            "coordinate:MI".parse::<SubsortOrder>()?,
            SubsortOrder::Coordinate(String::from("MI"))
        );

        assert_eq!(
            "unsorted:MI:coordinate".parse::<SubsortOrder>()?,
            SubsortOrder::Unsorted(String::from("MI:coordinate"))
        );

        assert!("".parse::<SubsortOrder>().is_err());
        assert!("noodles".parse::<SubsortOrder>().is_err());
        assert!("queryname".parse::<SubsortOrder>().is_err());
        assert!("QueryName".parse::<SubsortOrder>().is_err());

        Ok(())
    }
}

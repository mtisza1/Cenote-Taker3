Submit-block ::= {
  contact {
    contact {
      name name {
        last "Taker",
        first "Cenote",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "Cenote",
        div "Taker 2",
        city "Cenote",
        sub "DC",
        country "USA",
        street "1000 Circle St.",
        email "cenote@donotsubmit.net",
        postal-code "31111"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Taker",
            first "Cenote",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "Cenote",
        div "Taker 2",
        city "Cenote",
        sub "DC",
        country "USA",
        street "1000 Circle St.",
        postal-code "31111"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Taker",
              first "Cenote",
              middle "",
              initials "",
              suffix "",
              title ""
            }
          }
        }
      },
      title "Cenote-Taker 2 dummy template file. Do NOT submit"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "ALT EMAIL:cenote@donotsubmit.net"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title:None"
    }
  }
}
